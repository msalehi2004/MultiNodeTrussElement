#include "MultiNodeTruss2d.h"

#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

// Method to Read Command Arguments
void*
OPS_MultiNodeTruss2d()
{
	// check model dimension
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if (ndm != 2 || ndf != 3) {
		opserr << "WARNING! multiNodeTruss2d - ndm must be 2 and ndf must be 3\n";
		return 0;
	}

	// check number of arguments
	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	if (numRemainingArgs < 6) {
		opserr << "WARNING! multiNodeTruss2d - insufficient arguments\n";
		opserr << "Want: eleTag? nodeTags? matTag? A?\n";
		return 0;
	}

	// get inputs
	int eleTag;
	int numData = 1;
	if (OPS_GetIntInput(&numData, &eleTag) < 0) {
		opserr << "WARNING! multiNodeTruss2d - invalid element tag\n";
		return 0;
	}

	int subEleNo = numRemainingArgs - 3 - 1;
	numData = subEleNo + 1;
	int* nodeTags = new int[numData];
	if (OPS_GetIntInput(&numData, &nodeTags[0]) < 0) {
		opserr << "WARNING! multiNodeTruss2d - eleTag: " << eleTag << " - invalid node tags\n";
		return 0;
	}

	int matTag;
	numData = 1;
	if (OPS_GetIntInput(&numData, &matTag) < 0) {
		opserr << "WARNING! multiNodeTruss2d - eleTag: " << eleTag << " - invalid material tag\n";
		return 0;
	}

	UniaxialMaterial* axialMat = OPS_getUniaxialMaterial(matTag);
	if (!axialMat) {
		opserr << "WARNING! multiNodeTruss2d - eleTag: " << eleTag << " - material model not found\n";
		return 0;
	}

	double A;
	numData = 1;
	if (OPS_GetDoubleInput(&numData, &A) < 0) {
		opserr << "WARNING! multiNodeTruss2d - eleTag: " << eleTag << " - invalid A\n";
		return 0;
	}

	// create element
	Element* theEle = new MultiNodeTruss2d(eleTag, subEleNo, nodeTags, &axialMat, A);

	return theEle;
}

// Constructor 1 (for normal processing)
MultiNodeTruss2d::MultiNodeTruss2d(int tag, int N, int *nodes, UniaxialMaterial **mat, double a)
	: Element(tag, 0), connectedExternalNodes(N + 1),
	n(N), A(a), L(0.0), Lnew(0.0), K0(0), eps(0.0), F(0.0),
	theVectorSing(6), theVectorMult(3 * N + 3), theMatrixSing(6, 6), theMatrixMult(3 * N + 3, 3 * N + 3)

	// complete
{
	// Initialize theNodes Pointers
	theNodes = new Node *[n + 1];

	// Pointers to Nodes and Their IDs
	for (int i = 0; i < n + 1; i++) {
		connectedExternalNodes(i) = nodes[i];
		theNodes[i] = 0;
	}

	// Get Copy of Axial Material Model
	if (mat == 0) {
		opserr << "WARNING! MultiNodeTruss2d::MultiNodeTruss2d - element: " << this->getTag() << " - pointer to axial material model is null" << endln;
		exit(-1);
	}

	theMat = mat[0]->getCopy();

	if (!theMat) {
		opserr << "WARNING! MultiNodeTruss2d::MultiNodeTruss2d - element: " << this->getTag() << " - could not create copy of axial material model" << endln;
		exit(-1);
	}

	//this->revertToStart();
}

// Constructor 2 (for parallel processing)
MultiNodeTruss2d::MultiNodeTruss2d()
	: Element(0, 0), connectedExternalNodes(0),
	n(0), A(0.0), L(0.0), K0(0), eps(0.0), F(0.0)

	// complete
{

}

// Destructor
MultiNodeTruss2d::~MultiNodeTruss2d()
{
	// Delete theNodes
	if (theNodes)
		delete[] theNodes;

	// Delete Axial Material Pointer
	if (theMat)
		delete theMat;

	// Delete Initial Stiffness Pointer
	if (K0 != 0)
		delete K0;
}

// Definition of setDomain()
void
MultiNodeTruss2d::setDomain(Domain *theDomain)
{
	// Check Domain is not Null
	if (!theDomain) {
		for (int i = 0; i < n + 1; i++)
			theNodes[i] = 0;

		opserr << "ERROR! MultiNodeTruss2d::setDomain - element: " << this->getTag() << " - the domain is null\n";
		exit(0);
	}

	// Get Node Pointers
	for (int i = 0; i < n + 1; i++) {
		theNodes[i] = theDomain->getNode(connectedExternalNodes(i));

		if (!theNodes[i]) {
			opserr << "ERROR! MultiNodeTruss2d::setDomain - element: " << this->getTag() << " - node " << connectedExternalNodes(i) << " does not exist in the domain\n";
			exit(0);
		}

		// Check Node DOFs are 3
		if (theNodes[i]->getNumberDOF() != 3) {
			opserr << "ERROR! MultiNodeTruss2d::setDomain - element: " << this->getTag() << " - node " << connectedExternalNodes(i) << " has incorrect number of DOFs (not 3)\n";
			exit(0);
		}
	}

	// Get Initial Length
	if (this->initialLength() < 0)
		opserr << "ERROR! MultiNodeTruss2d::setDomain - element: " << this->getTag() << " - failed in initialLength\n";

	// Call DomainComponent Class Method
	this->DomainComponent::setDomain(theDomain);
}

// Definition of Methods Dealing with Nodes Information
int
MultiNodeTruss2d::getNumExternalNodes(void) const
{
	return (n + 1);
}

const ID &
MultiNodeTruss2d::getExternalNodes(void)
{
	return connectedExternalNodes;
}

Node **
MultiNodeTruss2d::getNodePtrs(void)
{
	return theNodes;
}

int
MultiNodeTruss2d::getNumDOF(void)
{
	return (3 * n + 3);
}

// Definition of Method Giving Lengths
int
MultiNodeTruss2d::initialLength()
{
	L = 0.0;

	for (int i = 0; i < n; i++) {
		const Vector &end1Crd = theNodes[i]->getCrds();
		const Vector &end2Crd = theNodes[i + 1]->getCrds();
		Vector dX = end2Crd - end1Crd;

		double L_i = dX.Norm();

		if (L_i < DBL_EPSILON) {
			opserr << "ERROR! MultiNodeTruss2d::initialLength - sub-element: " << i + 1 << " of element: " << this->getTag() << " - initial length is zero\n";
			return -1;
		}

		L += L_i;
	}

	return 0;
}

int
MultiNodeTruss2d::updateLength()
{
	Lnew = 0.0;

	for (int i = 0; i < n; i++) {
		const Vector &end1Crd = theNodes[i]->getCrds();
		const Vector &end2Crd = theNodes[i + 1]->getCrds();

		const Vector &dsp1 = theNodes[i]->getTrialDisp();
		const Vector &dsp2 = theNodes[i + 1]->getTrialDisp();

		Vector dX = end2Crd - end1Crd;

		dX(0) += (dsp2(0) - dsp1(0));
		dX(1) += (dsp2(1) - dsp1(1));

		double Ln_i = dX.Norm();

		if (Ln_i < DBL_EPSILON) {
			opserr << "ERROR! MultiNodeTruss2d::updateLength - sub-element: " << i + 1 << " of element: " << this->getTag() << " - new length is zero\n";
			return -1;
		}

		Lnew += Ln_i;
	}

	return 0;
}

// Definition of Methods Dealing with Transformation of Forces and Stiffnesses from Basic to Global System
Vector
MultiNodeTruss2d::forceTransf(int m)
{
	const Vector &end1Crd = theNodes[m]->getCrds();
	const Vector &end2Crd = theNodes[m + 1]->getCrds();

	const Vector &dsp1 = theNodes[m]->getTrialDisp();
	const Vector &dsp2 = theNodes[m + 1]->getTrialDisp();
	
	Vector dX = end2Crd - end1Crd;

	dX(0) += (dsp2(0) - dsp1(0));
	dX(1) += (dsp2(1) - dsp1(1));

	double Ln = dX.Norm();

	if (Ln < DBL_EPSILON) {
		opserr << "ERROR! MultiNodeTruss2d::initialTransf - sub-element: " << m << " of element: " << this->getTag() << " - new length is zero\n";
		exit(0);
	}

	const double cs = dX(0) / Ln;
	const double sn = dX(1) / Ln;

	theVectorSing(2) = theVectorSing(5) = 0.0;
	theVectorSing(0) = -cs;
	theVectorSing(1) = -sn;
	theVectorSing(3) = cs;
	theVectorSing(4) = sn;

	return theVectorSing;
}

Matrix
MultiNodeTruss2d::stiff(int m, double k, double f)
{
	const Vector &end1Crd = theNodes[m]->getCrds();
	const Vector &end2Crd = theNodes[m + 1]->getCrds();

	const Vector &dsp1 = theNodes[m]->getTrialDisp();
	const Vector &dsp2 = theNodes[m + 1]->getTrialDisp();

	Vector dX = end2Crd - end1Crd;

	dX(0) += (dsp2(0) - dsp1(0));
	dX(1) += (dsp2(1) - dsp1(1));

	double Ln = dX.Norm();

	if (Ln < DBL_EPSILON) {
		opserr << "ERROR! MultiNodeTruss2d::initialTransf - sub-element: " << m << " of element: " << this->getTag() << " - new length is zero\n";
		exit(0);
	}

	const double cs = dX(0) / Ln;
	const double sn = dX(1) / Ln;

	theMatrixSing.Zero();

	theMatrixSing(0, 0) = theMatrixSing(3, 3) = k * cs * cs;
	theMatrixSing(0, 1) = theMatrixSing(1, 0) = theMatrixSing(3, 4) = theMatrixSing(4, 3) = k * sn * cs;
	theMatrixSing(1, 1) = theMatrixSing(4, 4) = k * sn * sn;
	theMatrixSing(0, 3) = theMatrixSing(3, 0) = k * -(cs * cs);
	theMatrixSing(0, 4) = theMatrixSing(4, 0) = theMatrixSing(1, 3) = theMatrixSing(3, 1) = k * -(sn * cs);
	theMatrixSing(1, 4) = theMatrixSing(4, 1) = k * -(sn * sn);

	theMatrixSing(0, 0) += f / Ln;
	theMatrixSing(1, 1) += f / Ln;
	theMatrixSing(3, 3) += f / Ln;
	theMatrixSing(4, 4) += f / Ln;
	theMatrixSing(0, 3) += -f / Ln;
	theMatrixSing(3, 0) += -f / Ln;
	theMatrixSing(1, 4) += -f / Ln;
	theMatrixSing(4, 1) += -f / Ln;

	return theMatrixSing;
}

Matrix
MultiNodeTruss2d::initStiff(int m, double k)
{
	const Vector &end1Crd = theNodes[m]->getCrds();
	const Vector &end2Crd = theNodes[m + 1]->getCrds();

	Vector dX = end2Crd - end1Crd;

	double L0 = dX.Norm();

	if (L0 < DBL_EPSILON)
		opserr << "ERROR! MultiNodeTruss2d::initialTransf - sub-element: " << m << " of element: " << this->getTag() << " - new length is zero\n";

	const double cs = dX(0) / L0;
	const double sn = dX(1) / L0;

	theMatrixSing.Zero();

	theMatrixSing(0, 0) = theMatrixSing(3, 3) = k * cs * cs;
	theMatrixSing(0, 1) = theMatrixSing(1, 0) = theMatrixSing(3, 4) = theMatrixSing(4, 3) = k * sn * cs;
	theMatrixSing(1, 1) = theMatrixSing(4, 4) = k * sn * sn;
	theMatrixSing(0, 3) = theMatrixSing(3, 0) = k * -(cs * cs);
	theMatrixSing(0, 4) = theMatrixSing(4, 0) = theMatrixSing(1, 3) = theMatrixSing(3, 1) = k * -(sn * cs);
	theMatrixSing(1, 4) = theMatrixSing(4, 1) = k * -(sn * sn);

	return theMatrixSing;
}

// Definition of Methods Dealing with Element's State
int
MultiNodeTruss2d::commitState(void)
{
	int err = 0;

	// Element commitState()
	if ((err = this->Element::commitState())) {
		opserr << "ERROR! MultiNodeTruss2d::commitState - element: " << this->getTag() << " - failed in base class\n";
	}

	// Commit Material State Variables
	err += theMat->commitState();

	// complete committing the variables

	return err;
}

int
MultiNodeTruss2d::revertToLastCommit(void)
{
	int err = 0;

	// Revert Material State Variables to Last Committed State
	err += theMat->revertToLastCommit();

	// complete reverting the variables

	return err;
}

int
MultiNodeTruss2d::revertToStart(void)
{
	int err = 0;

	// Revert Material State Variables to Start
	err += theMat->revertToStart();

	// complete committing the variables

	return err;
}

// Definition of Methods Dealing with Matrix/Vector Assembling
void
MultiNodeTruss2d::add2Matrix(Matrix &A, const Matrix &B, int rowStart, int rowEnd, int colStart, int colEnd, double fact)
{
	int rowsNo = rowEnd - rowStart + 1;
	int colsNo = colEnd - colStart + 1;

	if (B.noRows() != rowsNo)
		opserr << "ERROR! MultiNodeTruss2d::assembleMatrix - element: " << this->getTag() << " - incompatible number of rows to assemble\n";

	if (B.noCols() != colsNo)
		opserr << "ERROR! MultiNodeTruss2d::assembleMatrix - element: " << this->getTag() << " - incompatible number of columns to assemble\n";

	if ((A.noRows() - 1) < rowEnd)
		opserr << "ERROR! MultiNodeTruss2d::assembleMatrix - element: " << this->getTag() << " - receiving matrix has less rows than needed\n";

	if ((A.noCols() - 1) < colEnd)
		opserr << "ERROR! MultiNodeTruss2d::assembleMatrix - element: " << this->getTag() << " - receiving matrix has less columns than needed\n";

	int i = 0;
	int j = 0;

	for (int row = rowStart; row <= rowEnd; row++) {
		for (int col = colStart; col <= colEnd; col++) {
			A(row, col) += (fact * B(i, j));
			j++;
		}

		j = 0;
		i++;
	}
}

void
MultiNodeTruss2d::add2Vector(Vector &A, const Vector &B, int rowStart, int rowEnd, double fact)
{
	int rowsNo = rowEnd - rowStart + 1;

	if (B.Size() != rowsNo)
		opserr << "ERROR! MultiNodeTruss2d::assembleVector - element: " << this->getTag() << " - incompatible number of rows to assemble\n";

	if ((A.Size() - 1) < rowEnd)
		opserr << "ERROR! MultiNodeTruss2d::assembleVector - element: " << this->getTag() << " - receiving matrix has less rows than needed\n";

	int i = 0;

	for (int row = rowStart; row <= rowEnd; row++) {
		A(row) += (fact * B(i));

		i++;
	}
}

// Element Solution
int
MultiNodeTruss2d::update(void)
{
	// Compute Current Strain
	if (this->updateLength() < 0) {
		opserr << "ERROR! MultiNodeTruss2d::update - element: " << this->getTag() << " - failed to update length\n";
		return -1;
	}
	eps = Lnew / L - 1.0;

	// Compute Current Axial Force
	theMat->setTrialStrain(eps);
	F = A * (theMat->getStress());

	return 0;
}

// Definition of Methods Used to Determine Stiffness/Mass Matrices
const Matrix &
MultiNodeTruss2d::getTangentStiff(void)
{
	// Determine Element Axial Stiffness in Basic System
	double k = A * (theMat->getTangent()) / L;

	// Get Tranformed Stiffnesses for Subelements and Assemble Entire Stiffness Matrix
	theMatrixMult.Zero();

	for (int i = 0; i < n; i++)
		this->add2Matrix(theMatrixMult, this->stiff(i, k, F), 3 * i, 3 * i + 5, 3 * i, 3 * i + 5, 1.0);

	return theMatrixMult;
}

const Matrix &
MultiNodeTruss2d::getInitialStiff(void)
{
	// Check for Quick Return
	if (K0 != 0)
		return *K0;

	// Determine Element's Initial Axial Stiffness in Basic System
	double k0 = A * (theMat->getInitialTangent()) / L;

	// Get Tranformed Stiffnesses for Subelements and Assemble Entire Stiffness Matrix
	theMatrixMult.Zero();

	for (int i = 0; i < n; i++)
		this->add2Matrix(theMatrixMult, this->initStiff(i, k0), 3 * i, 3 * i + 5, 3 * i, 3 * i + 5, 1.0);

	K0 = new Matrix(theMatrixMult);

	return *K0;
}

const Matrix &
MultiNodeTruss2d::getMass(void)
{
	theMatrixMult.Zero();

	return theMatrixMult;
}

// Definition of Methods to Get Forces
const Vector &
MultiNodeTruss2d::getResistingForce(void)
{
	// Get Tranformed Forces for Subelements and Assemble Entire Force Vector
	theVectorMult.Zero();

	for (int i = 0; i < n; i++)
		this->add2Vector(theVectorMult, this->forceTransf(i), 3 * i, 3 * i + 5, F);

	return theVectorMult;
}

const Vector &
MultiNodeTruss2d::getResistingForceIncInertia()
{
	// Compute the current resisting force
	theVectorMult = this->getResistingForce();

	// add the damping forces if rayleigh damping
	if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		theVectorMult += this->getRayleighDampingForces();

	return theVectorMult;
}

// Definition of Print Command
void
MultiNodeTruss2d::Print(OPS_Stream &s, int flag)
{
	s << "Element Tag: " << this->getTag() << endln;
	s << "Type: MultiNodeTruss2d" << endln;
	s << "Connected Node Tags: firstNode " << connectedExternalNodes(0)
		<< ", lastNode " << connectedExternalNodes(n) << endln;
	s << "Cross Section Area: " << A << endln;
	s << "Axial Material Model Tag: " << theMat->getTag() << endln;
}

// Definition of Response Parameters
Response *
MultiNodeTruss2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	// Define and Initialize theResponse
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", this->getClassType());
	output.attr("eleTag", this->getTag());
	output.attr("nodeFirst", connectedExternalNodes[0]);
	output.attr("nodeLast", connectedExternalNodes[n]);

	// Global Forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, theVectorMult);
	}

	// Basic Forces
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "axialForce") == 0) {
		output.tag("ResponseType", "N_J");

		theResponse = new ElementResponse(this, 2, F);
	}

	// Axial Deformation
	else if (strcmp(argv[0], "deformation") == 0 || strcmp(argv[0], "axialDisplacement") == 0) {
		output.tag("ResponseType", "ux_J");

		theResponse = new ElementResponse(this, 3, eps * L);
	}

	// Axial Strain
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "axialStrain") == 0) {
		output.tag("ResponseType", "eps");

		theResponse = new ElementResponse(this, 4, eps);
	}

	// Slackness
	//else if (strcmp(argv[0], "slackness") == 0 || strcmp(argv[0], "slackLength") == 0) {
		//output.tag("ResponseType", "ux_J");

		//theResponse = new ElementResponse(this, 5, theMat->getSlackStrain() * L);
	//}

	// a material quantity
	else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "-material") == 0) {

		theResponse = theMat->setResponse(&argv[1], argc - 1, output);
	}

	return theResponse;
}

// Definition of Method to Get Response Parameters
int
MultiNodeTruss2d::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1: // Global Forces
		return eleInfo.setVector(this->getResistingForce());

	case 2: // Basic Forces
		return eleInfo.setDouble(F);

	case 3: // Axial Deformation
		return eleInfo.setDouble(eps * L);

	case 4: // Axial Strain
		return eleInfo.setDouble(eps);

	//case 5: // Slackness
		//return eleInfo.setDouble(theMat->getSlackStrain() * L);

	default:
		return -1;
	}
}

// Definition of Method to Display Element
int
MultiNodeTruss2d::displaySelf(Renderer &theViewer, int displayMode, float fact,
	const char** displayModes, int numModes)
{
	int err = 0;

	static Vector v1(3);
	static Vector v2(3);

	for (int ele = 0; ele < n; ele++) {
		// first determine the end points of the beam based on
		// the display factor (a measure of the distorted image)

		if (displayMode >= 0) {
			theNodes[ele]->getDisplayCrds(v1, fact);
			theNodes[ele + 1]->getDisplayCrds(v2, fact);

			err += theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
		}
		else {
			theNodes[ele]->getDisplayCrds(v1, 0.0);
			theNodes[ele + 1]->getDisplayCrds(v2, 0.0);

			int mode = displayMode  *  -1;

			const Matrix &eigen1 = theNodes[ele]->getEigenvectors();
			const Matrix &eigen2 = theNodes[ele + 1]->getEigenvectors();

			if (eigen1.noCols() >= mode) {
				for (int i = 0; i < 2; i++) {
					v1(i) += eigen1(i, mode - 1)*fact;
					v2(i) += eigen2(i, mode - 1)*fact;
				}
			}

			err += theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
		}
	}

	return err;
}

// Definition of Methods Dealing with Parallel Processing
int
MultiNodeTruss2d::sendSelf(int commitTag, Channel &theChannel)
{
	opserr << "ERROR! MultiNodeTruss2d::sendSelf - element: " << this->getTag() << " - incapable of parallel processing\n";
	return -1;
}

int
MultiNodeTruss2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	opserr << "ERROR! MultiNodeTruss2d::recvSelf - element: " << this->getTag() << " - incapable of parallel processing\n";
	return -1;
}
