#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction
class FENewActiveFiberContraction : public FEMaterial
{
public:
	FENewActiveFiberContraction(FEModel* pfem);

	//! initialization
	void Init();

	//! calculate the fiber stress
	double FiberStress(double lamd);

	//! active contraction stiffness contribution
	double FiberStiffness(double lamd);

public:
	//! set the property attributes
	bool SetAttribute(const char* szname, const char* szval);

	//! Get the activation
	double GetActivation() { return m_ascl; }

protected:
	double	m_ascl;		//!< activation scale factor
	double	m_smax;
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
// 	double	m_Tmax;		//!< activation scale factor
// 	double	m_ca0;		//!< intracellular calcium concentration
// 	double	m_camax;	//!< peak calcium concentration
// 	double	m_beta;		//!< shape of peak isometric tension-sarcomere length relation
// 	double	m_l0;		//!< unloaded length
// 	double	m_refl;		//!< sarcomere length

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FENewFiberMaterial : public FEMaterial
{
public:
	//! Constructor
	FENewFiberMaterial(FEModel* pfem);

	//! Initialization
	void Init();

	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);

public:
	void Serialize(DumpFile& ar);

	//! Set the active contraction property
	void SetActiveContraction(FENewActiveFiberContraction* pma) { m_pafc = pma; }

	//! get the active contraction property
	FEMaterial* GetActiveContraction() { return m_pafc; }

	//! get activation
	double GetActivation() { return (m_pafc ? m_pafc->GetActivation() : 0.0); }

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers

	double	m_lam1;		//!< fiber stretch for straightened fibers

	//--- time varying elastance active contraction data ---
	FENewActiveFiberContraction*	m_pafc;

	DECLARE_PARAMETER_LIST();
};
