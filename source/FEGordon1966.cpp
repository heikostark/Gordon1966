#include "stdafx.h"
#include "FEGordon1966.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEGordon1966, FEUncoupledMaterial)
	ADD_PARAMETER(c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(c2, FE_PARAM_DOUBLE, "c2");
	
	ADD_PARAMETER(m_fib.m_c3, FE_PARAM_DOUBLE, "c3");
	ADD_PARAMETER(m_fib.m_c4, FE_PARAM_DOUBLE, "c4");
	ADD_PARAMETER(m_fib.m_c5, FE_PARAM_DOUBLE, "c5");
	ADD_PARAMETER(m_fib.m_lam1, FE_PARAM_DOUBLE, "lam_max");
	
	ADD_PARAMETER(m_fib.m_pafc->m_ascl , FE_PARAM_DOUBLE, "ascl" );	
	ADD_PARAMETER(m_fib.m_pafc->m_smax , FE_PARAM_DOUBLE, "smax" );
	ADD_PARAMETER(m_fib.m_pafc->m_ax , FE_PARAM_DOUBLE, "ax" );
	ADD_PARAMETER(m_fib.m_pafc->m_ay , FE_PARAM_DOUBLE, "ay" );
	ADD_PARAMETER(m_fib.m_pafc->m_bx , FE_PARAM_DOUBLE, "bx" );
	ADD_PARAMETER(m_fib.m_pafc->m_by , FE_PARAM_DOUBLE, "by" );
	ADD_PARAMETER(m_fib.m_pafc->m_cx , FE_PARAM_DOUBLE, "cx" );
	ADD_PARAMETER(m_fib.m_pafc->m_cy , FE_PARAM_DOUBLE, "cy" );
	ADD_PARAMETER(m_fib.m_pafc->m_dx , FE_PARAM_DOUBLE, "dx" );
	ADD_PARAMETER(m_fib.m_pafc->m_dy , FE_PARAM_DOUBLE, "dy" );
	ADD_PARAMETER(m_fib.m_pafc->m_ex , FE_PARAM_DOUBLE, "ex" );
	ADD_PARAMETER(m_fib.m_pafc->m_ey , FE_PARAM_DOUBLE, "ey" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEGordon1966::FEGordon1966(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem)  {   printf("Gordon1966\n"); }

//-----------------------------------------------------------------------------
// Data initialization
void FEGordon1966::Init()
{
	m_fib.Init();	// first init for m_fib.m_pafc
	FEUncoupledMaterial::Init();
}

//-----------------------------------------------------------------------------
// This material has two properties (the fiber material and the active contraction material)
int FEGordon1966::Properties() { return 2; }

//-----------------------------------------------------------------------------
FECoreBase* FEGordon1966::GetProperty(int n)
{
	if (n == 0) return &m_fib;
	if (n == 1) return m_fib.GetActiveContraction();
	return 0;
}

//-----------------------------------------------------------------------------
//! find a material property index ( returns <0 for error)
int FEGordon1966::FindPropertyIndex(const char* szname)
{
	return 1;
}

//-----------------------------------------------------------------------------
//! set a material property (returns false on error)
bool FEGordon1966::SetProperty(int i, FECoreBase* pm)
{
	if (i == 1)
	{
		FENewActiveFiberContraction* pma = dynamic_cast<FENewActiveFiberContraction*>(pm);
		if (pma) { m_fib.SetActiveContraction(pma); return true; }
	}
	return false;
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file 
void FEGordon1966::Serialize(DumpFile &ar)
{
	// serialize the base class parameters
	FEUncoupledMaterial::Serialize(ar);
	// serialize fiber data
	m_fib.Serialize(ar);
}

//-----------------------------------------------------------------------------
mat3ds FEGordon1966::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = c1;
	double W2 = c2;
	// ------------------------------------------------

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	// calculate stress s = pI + 2/J * dev(T) 
	mat3ds s = T.dev()*(2.0/J);

	// add the fiber stress
	s += m_fib.Stress(mp);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEGordon1966::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0/J;

	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B*B;

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2;
	W1 = c1;
	W2 = c2;
	// ------------------------------------

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// mean pressure
	//double p = pt.avgp;

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);
	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c + m_fib.Tangent(mp);
}
