#pragma once
#include "FEBioMech/FEUncoupledMaterial.h"
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

	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2
	
	FENewFiberMaterial	m_fib;
	
	//! Initialization
	bool Init();

	//! serialize material data
	void Serialize(DumpStream& ar);
	
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	// declare parameter list
	DECLARE_PARAMETER_LIST(); // virtual void BuildParamList();
};
