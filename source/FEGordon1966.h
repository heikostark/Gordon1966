#pragma once
#include "FEBioMech/FEUncoupledMaterial.h"
//#include "FEBioMech/FEUncoupledFiberExpLinear.h"
//#include "FEBioMech/FEActiveContractionMaterial.h"
#include <FECore/FEModelParam.h>

#include "FENewFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Transversely Isotropic Mooney-Rivlin material

//! This material has an isotopric Mooney-Rivlin basis and single preferred
//! fiber direction.

class FEGordon1966: public FEUncoupledMaterial
{
public:
	//! constructor
	FEGordon1966(FEModel* pfem);

private:
	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2
	
	FENewFiberMaterial	m_fib;
	
public:	
	//! Initialization
	virtual bool Init();

	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

    //! create material point data
    FEMaterialPointData* CreateMaterialPointData() override;
	
protected:
	//FEFiberExpLinearUC		       m_fib;
    //FEActiveContractionMaterial*   m_ac;
	FEVec3dValuator* m_fiber;
	
	// This macro defines that the class will define a material parameter list.
	// The material parameter list itself is defined elsewhere (e.g. in the .cpp file.)
	DECLARE_FECORE_CLASS();
};
