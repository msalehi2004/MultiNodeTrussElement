#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include "MultiNodeTruss2d.h"
#include "MultiNodeTruss3d.h"

extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addMultiNodeTruss(ClientData clientData,
	Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
	TclModelBuilder *theTclBuilder)
{
	int eleArgStart = 1;

	// ensure the destructor has not been called
	if (theTclBuilder == 0)  {
		opserr << "WARNING builder has been destroyed - MultiNodeTruss\n";
		return TCL_ERROR;
	}

	Element *theElement = 0;
	int ndm = theTclBuilder->getNDM();
	int ndf = theTclBuilder->getNDF();
	int tag;

	if (ndm == 2)  {
		// check plane frame problem has 3 dof per node
		if (ndf != 3)  {
			opserr << "WARNING invalid ndf: " << ndf;
			opserr << ", for plane problem need 3 - MultiNodeTruss\n";
			return TCL_ERROR;
		}

		// check the number of arguments is correct
		if ((argc - eleArgStart) < 6)  {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: multiNodeTruss eleTag? nodeTags? matTag? A?\n";
			return TCL_ERROR;
		}

		// get inputs
		int eleNo = argc - eleArgStart - 4 - 1;
		int matTag;
		double A;

		if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK)  {
			opserr << "WARNING invalid MultiNodeTruss eleTag";
			return TCL_ERROR;
		}

		int *nodeTags = new int[eleNo + 1];

		for (int i = 0; i < (eleNo + 1); i++) {
			if (Tcl_GetInt(interp, argv[2 + eleArgStart + i], &nodeTags[i]) != TCL_OK)  {
				opserr << "WARNING invalid " << i+1 << "-th Node";
				opserr << " - MultiNodeTruss element: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetInt(interp, argv[2 + eleArgStart + eleNo + 1], &matTag) != TCL_OK)  {
			opserr << "WARNING invalid matTag";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		UniaxialMaterial* axialMat = OPS_getUniaxialMaterial(matTag);

		if (!axialMat)  {
			opserr << "WARNING material model not found";
			opserr << " - uniaxialMaterial: " << matTag;
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3 + eleArgStart + eleNo + 1], &A) != TCL_OK)  {
			opserr << "WARNING invalid A";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		// now create the MultiNodeTruss
		theElement = new MultiNodeTruss2d(tag, eleNo, nodeTags, &axialMat, A);

		if (!theElement)  {
			opserr << "WARNING ran out of memory creating element";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		// then add the MultiNodeTruss to the domain
		if (theTclDomain->addElement(theElement) == false)  {
			opserr << "WARNING could not add element to the domain";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			delete theElement;
			return TCL_ERROR;
		}
	}
	else if (ndm == 3)  {
		// check plane frame problem has 6 dof per node
		if (ndf != 6)  {
			opserr << "WARNING invalid ndf: " << ndf;
			opserr << ", for plane problem need 6 - MultiNodeTruss\n";
			return TCL_ERROR;
		}

		// check the number of arguments is correct
		if ((argc - eleArgStart) < 6)  {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: multiNodeTruss eleTag? nodeTags? matTag? A?\n";
			return TCL_ERROR;
		}

		// get the id and end nodes
		int eleNo = argc - eleArgStart - 4 - 1;
		int matTag;
		double A;

		if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK)  {
			opserr << "WARNING invalid MultiNodeTruss eleTag";
			return TCL_ERROR;
		}

		int *nodeTags = new int[eleNo + 1];

		for (int i = 0; i < (eleNo + 1); i++) {
			if (Tcl_GetInt(interp, argv[2 + eleArgStart + i], &nodeTags[i]) != TCL_OK)  {
				opserr << "WARNING invalid " << i + 1 << "-th Node";
				opserr << " - MultiNodeTruss element: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetInt(interp, argv[2 + eleArgStart + eleNo + 1], &matTag) != TCL_OK)  {
			opserr << "WARNING invalid matTag";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		UniaxialMaterial* axialMat = OPS_getUniaxialMaterial(matTag);

		if (!axialMat)  {
			opserr << "WARNING material model not found";
			opserr << " - uniaxialMaterial: " << matTag;
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3 + eleArgStart + eleNo + 1], &A) != TCL_OK)  {
			opserr << "WARNING invalid A";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		// now create the MultiNodeTruss
		theElement = new MultiNodeTruss3d(tag, eleNo, nodeTags, &axialMat, A);

		if (!theElement)  {
			opserr << "WARNING ran out of memory creating element";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			return TCL_ERROR;
		}

		// then add the MultiNodeTruss to the domain
		if (theTclDomain->addElement(theElement) == false)  {
			opserr << "WARNING could not add element to the domain";
			opserr << " - MultiNodeTruss element: " << tag << endln;
			delete theElement;
			return TCL_ERROR;
		}
	}
	else  {
		opserr << "WARNING multiNodeTruss command only works when ndm is 2 or 3, ndm: ";
		opserr << ndm << endln;
		return TCL_ERROR;
	}
	
	// if get here we have sucessfully created the HSRjoint and added it to the domain
	return TCL_OK;
}

