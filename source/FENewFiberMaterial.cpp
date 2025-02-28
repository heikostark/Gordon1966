#include "stdafx.h"
#include "FENewFiberMaterial.h"
#include "FEBioMech/FEElasticMaterial.h"
#include <stdlib.h>
#ifdef HAVE_GSL
#include "gsl/gsl_sf_expint.h"
#endif

//-----------------------------------------------------------------------------
FENewActiveFiberContraction::FENewActiveFiberContraction(FEModel* pfem) : FEMaterial(pfem)
{
	m_ascl = 0;
	m_smax = 1;
	m_ax = 0.357; 
	m_ay = 0;
	m_bx = 0.596;
	m_by = 0.84;
	m_cx = 0.929;
	m_cy = 1;
	m_dx = 1;
	m_dy = 1;
	m_ex = 1.521;
	m_ey = 0;
	printf("    NewActiveFiberContraction\n");
}

//-----------------------------------------------------------------------------
void FENewActiveFiberContraction::Init() { }

//-----------------------------------------------------------------------------
double FENewActiveFiberContraction::FiberStress(double lamd)
{
	// --- active contraction contribution ---
	double saf = 0.0;
	if (m_ascl > 0)
	{
		double ctenslm = m_ascl; // activation scale factor
		double Wl = 0;
		/*if (lamd < m_ax) Wl = 0;
		else*/ if (lamd < m_bx) Wl = m_ay + (lamd-m_ax)*((m_by-m_ay)/(m_bx-m_ax));
		else if (lamd < m_cx) Wl = m_by + (lamd-m_bx)*((m_cy-m_by)/(m_cx-m_bx));
		else if (lamd < m_dx) Wl = m_cy + (lamd-m_cx)*((m_dy-m_cy)/(m_dx-m_cx));
		else if (lamd < m_ex) Wl = m_dy + (lamd-m_dx)*((m_ey-m_dy)/(m_ex-m_dx));
		/*else Wl = 0;*/

		saf = ctenslm * (m_smax * Wl); // activation * (max stress * W on lamd)
	}
	return saf;
}

//-----------------------------------------------------------------------------
double FENewActiveFiberContraction::FiberStiffness(double lamd)
{
/*	if (lcna >= 0)
	{
		double ctenslm = m_plc->Value();

		// current sarcomere length
		double strl = refl*lamd;

		// sarcomere length change
		double dl = strl - l0;

		if (dl >= 0) W44 += J*2*beta*refl*exp(-beta*dl);
	}
*/
	return 0.0;
}

//=============================================================================
//=============================================================================
//=============================================================================

//-----------------------------------------------------------------------------
FENewFiberMaterial::FENewFiberMaterial(FEModel* pfem) : FEMaterial(pfem)
{
	m_c3 = m_c4 = m_c5 = 0;
	m_lam1 = 1;
	SetActiveContraction (new FENewActiveFiberContraction(pfem));
	printf("  NewFiberMaterial\n");	
}

//-----------------------------------------------------------------------------
void FENewFiberMaterial::Init()
{        
	if (m_pafc) m_pafc->Init();
}

//-----------------------------------------------------------------------------
// Fiber material stress
mat3ds FENewFiberMaterial::Stress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0 / J;
	double Jm13 = pow(J, -1.0/3.0);
	double twoJi = 2.0*Ji;

	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = 0;
	if (lamd > 1)
	{
		double lamdi = 1.0/lamd;
		double Wl;
		if (lamd < m_lam1)
		{
			Wl = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1-1))-1) - m_c5*m_lam1;
			Wl = lamdi*(m_c5*lamd + c6);
		}
		W4  = 0.5*lamdi*Wl;
	}
	else 
	{
		W4 = 0;
	}	

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// ---
	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);
	
	// calculate stress: 
	mat3ds s = T.dev()*twoJi;

	// --- active contraction contribution ---
	if (m_pafc) s += AxA*m_pafc->FiberStress(lamd);

	return s;
}

//-----------------------------------------------------------------------------
// Fiber material tangent
tens4ds FENewFiberMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
	double Jm23 = Jm13*Jm13;
	double Ji = 1.0/J;

	// get initial local material axis
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	double W4, W44;
	if (lamd >= 1)
	{
		double lamdi = 1.0/lamd;
		double Wl, Wll;
		if (lamd < m_lam1)
		{
			Wl  = lamdi*m_c3*(exp(m_c4*(lamd - 1)) - 1);
			Wll = m_c3*lamdi*(m_c4*exp(m_c4*(lamd - 1)) - lamdi*(exp(m_c4*(lamd-1))-1));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1-1))-1) - m_c5*m_lam1;
			Wl  = lamdi*(m_c5*lamd + c6);
			Wll = -c6*lamdi*lamdi;
		}
		W4  = 0.5*lamdi*Wl;
		W44 = 0.25*lamdi*lamdi*(Wll - lamdi*Wl);
	}
	else 
	{
		W4 = 0;
		W44 = 0;
	}

	// --- add active contraction stiffness ---
	if (m_pafc) W44 += m_pafc->FiberStiffness(lamd);

	// --- calculate tangent ---

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4  = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0/3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI/3.0)*(4.0/3.0*Ji*WC) + IxI*(4.0/9.0*Ji*CWWC) + cw;

	return c;
}

//-----------------------------------------------------------------------------
// Fiber material strain energy density
//
double FENewFiberMaterial::StrainEnergyDensity(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
    
	// get the initial fiber direction
	vec3d a0;
	a0.x = pt.m_Q[0][0];
	a0.y = pt.m_Q[1][0];
	a0.z = pt.m_Q[2][0];
    
	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;
    
	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde
    
	// strain energy density
	double sed = 0.0;
#ifdef HAVE_GSL
	if (lamd > 1)
	{
		if (lamd < m_lam1)
		{
			sed = m_c3*(exp(-m_c4)*
                        (gsl_sf_expint_Ei(m_c4*lamd)-gsl_sf_expint_Ei(m_c4))
                        - log(lamd));
		}
		else
		{
			double c6 = m_c3*(exp(m_c4*(m_lam1-1))-1) - m_c5*m_lam1;
			sed = m_c5*(lamd-1) + c6*log(lamd);
		}
	}
#endif
	// --- active contraction contribution to sed is zero ---
    
	return sed;
}

//-----------------------------------------------------------------------------
void FENewFiberMaterial::Serialize(DumpFile& ar)
{
	FEMaterial::Serialize(ar);
	if (ar.IsSaving())
	{
		if (m_pafc)
		{
			ar << 1;
			m_pafc->Serialize(ar);
		}
		else ar << 0;
	}
	else
	{
		int nafc;
		ar >> nafc;
		if (nafc == 1)
		{
			m_pafc = new FENewActiveFiberContraction(GetFEModel());
			m_pafc->Serialize(ar);
		}
	}
}