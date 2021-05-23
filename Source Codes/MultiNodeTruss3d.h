#ifndef MULTINODETRUSS3D_H
#define MULTINODETRUSS3D_H

// Directives
#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <UniaxialMaterial.h>

// Forward Declaration
class Response;

class MultiNodeTruss3d : public Element
{
public:
	// Constructors
	MultiNodeTruss3d();
	MultiNodeTruss3d(int tag, int N,  int *nodes, UniaxialMaterial **mat, double a);

	// Destructor
	~MultiNodeTruss3d();

	// Method to Get Class Type
	const char *getClassType() const { return "MultiNodeTruss3d"; };

	// Method to Initialize the Domain; base class: DomainComponent
	void setDomain(Domain *theDomain);

	// Methods to Obtain Information about DOFs and Connectivity; base class: Element
	int getNumExternalNodes(void) const;
	const ID &getExternalNodes(void);
	Node **getNodePtrs(void);
	int getNumDOF(void);

	// Methods to Set the State of Element; base class: Element
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	int update(void);

	// Methods to Obtain Stiffness/Mass Matrices; base class: Element
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getMass(void);

	// Methods to Obtain Resisting Forces; base class: Element
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// Methods to Obtain Information Specific to Element; base class: Element
	void Print(OPS_Stream &s, int flag = 0);
	Response *setResponse(const char **argv, int argc, OPS_Stream &output);
	int getResponse(int responseID, Information &eleInfo);

	// Method to Display Element
	int	displaySelf(Renderer &theViewer, int displayMode, float fact, const char** displayModes = 0, int numModes = 0);

	// Methods to Do Parallel Processing; base class: Channel
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

private:
	// Private Attributes
	ID connectedExternalNodes;              // contains tags of connecting nodes
	Node **theNodes;                        // pointer to nodes
	UniaxialMaterial *theMat;				// pointer to uniaxial material

	int n;									// number of elements

	double L;								// element's initial length
	double A;								// cross section area

	Matrix glTrans;							// global system to local system transformation matrix
	Matrix lbTrans;						    // local system to basic system transformation matrix
	Matrix *K0;								// pointer to initial stiffness matrix in the basic system

	// State Variables
	double F;								// constant axial force in the truss elements
	double eps;								// constant axial strain in the truss elements

	// complete the state variables

	// Private Methods
	void add2Matrix(Matrix &A, const Matrix &B, int rowStart, int rowEnd, int colStart, int colEnd, double fact);
	void add2Vector(Vector &A, const Vector &B, int rowStart, int rowEnd, double fact);

	double initialLength(void);
	double newLength(void);

	Vector forceTransf(int m);
	Matrix stiff(int m, double k, double f);
	Matrix initStiff(int m, double k);

	// Static Class-Wide Variables
	Matrix theMatrixSing;
	Matrix theMatrixMult;
	Vector theVectorSing;
	Vector theVectorMult;
};

#endif // MultiNodeTruss3d_H
