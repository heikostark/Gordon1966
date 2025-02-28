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

	double	c1;	//!< Mooney-Rivlin coefficient C1
	double	c2;	//!< Mooney-Rivlin coefficient C2
	
	FENewFiberMaterial	m_fib;
	
	//! Initialization
	void Init();

	//! serialize material data
	void Serialize(DumpFile& ar);

	//! return the number of properties
	int Properties();

	//! return a pointer to the property
	FECoreBase* GetProperty(int n);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FECoreBase* pm);

	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	// declare parameter list
	DECLARE_PARAMETER_LIST(); // virtual void BuildParamList();
};
