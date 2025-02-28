#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction
class FENewActiveFiberContraction : public FEMaterial
{
public:
	FENewActiveFiberContraction(FEModel* pfem);

	double	m_ascl;		//!< activation scale factor - lc=".."
	double	m_smax;		//!< max activation
//! force-length curve	
	double	m_ax;
	double	m_ay;
	double	m_bx;
	double	m_by;
	double	m_cx;
	double	m_cy;
	double	m_dx;
	double	m_dy;
	double	m_ex;
	double	m_ey;	
	
	//! initialization
	void Init();

	//! calculate the fiber stress
	double FiberStress(double lamd);

	//! active contraction stiffness contribution
	double FiberStiffness(double lamd);

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FENewFiberMaterial : public FEMaterial
{
public:
	//! Constructor
	FENewFiberMaterial(FEModel* pfem);

	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers
	double	m_lam1;		//!< fiber stretch for straightened fibers

	//--- time varying elastance active contraction data ---
	FENewActiveFiberContraction*	m_pafc; // pointer to contraction
	
	//! Initialization
	void Init();

	//! Set the active contraction property
	void SetActiveContraction(FENewActiveFiberContraction* pma) { m_pafc = pma; }

	//! get the active contraction property
	FEMaterial* GetActiveContraction() { return m_pafc; }

	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	void Serialize(DumpFile& ar);
	
	DECLARE_PARAMETER_LIST();
};