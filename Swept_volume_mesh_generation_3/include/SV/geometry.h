#ifndef SV_GEOMETRY_H 
#define SV_GEOMETRY_H

#include <geometry/triangle.h>
#include <math/vector3.h>

namespace SV{

	
bool inline teIntersect(const leprechaun::geometry::Triangle<float>& T, const Vector3<float>& p, const Vector3<float>& q)
{
	// Schnitttest auf Grundlage des Separating Axis Thm
	// mit early out
	const leprechaun::geometry::Triangle<float>& T1 = T;
	Vector3<float> n,u,v;
	float d,dmin,dmax,emin,emax;

	// Achse 1					
	n = T1.normal();  					// Richtung der ersten Dreiecksnormalen
	dmin = dmax = n*(p-T1[0]);
	d = n*(q-T1[0]);
	dmin = std::min(dmin,d);
	dmax = std::max(dmax,d);

	if (dmin > 0) return false;
	if (dmax < 0) return false;

	// Achsen 3-11	
	for(int i=0;i<3;i++) 
	{
		u = T1[i]-T1[(i+1)%3];
		v = q - p;
		n = u ^ v;
		dmin = dmax = 0;
		for(int k=1;k<3;k++) {
			d = n*(T1[k]-T1[0]);
			dmin = std::min(dmin,d);
			dmax = std::max(dmax,d);
		}
		emin = emax = n*(p-T1[0]);
		d = n*(q-T1[0]);
		emin = std::min(emin,d);
		emax = std::max(emax,d);
		if (emin > dmax) return false;
		if (emax < dmin) return false;
	}
	return true;
}

inline float vvDist2( const Vector3<float>& u, const Vector3<float>& v )
{
	  return (u[0]-v[0])*(u[0]-v[0])+(u[1]-v[1])*(u[1]-v[1])+(u[2]-v[2])*(u[2]-v[2]);
}

//float veDist2(const float p[], const float u[], const float v[]) 
inline float veDist2( const Vector3<float>& p, const Vector3<float>& u, const Vector3<float>& v )
{
	Vector3<float> pu,vu;
	float puvu,vuvu,pupu;

	pu[0] = p[0]-u[0]; pu[1] = p[1]-u[1]; pu[2] = p[2]-u[2];
	vu[0] = v[0]-u[0]; vu[1] = v[1]-u[1]; vu[2] = v[2]-u[2];
	vuvu = vu[0]*vu[0]+vu[1]*vu[1]+vu[2]*vu[2];
	puvu = pu[0]*vu[0]+pu[1]*vu[1]+pu[2]*vu[2];
	pupu = pu[0]*pu[0]+pu[1]*pu[1]+pu[2]*pu[2];

	if ((puvu > 0.0) & (puvu < vuvu)) return pupu-puvu*puvu/vuvu;
	return 1e10;
}
float eeDist2( const Vector3<float>& a, const Vector3<float>& b, const Vector3<float>& u, const Vector3<float>& v )
{
	Vector3<float>  ua,vu,ba,n,nxua;
	float nn,nua,l1,l2;

	ba[0] = b[0]-a[0]; ba[1] = b[1]-a[1]; ba[2] = b[2]-a[2];
	ua[0] = u[0]-a[0]; ua[1] = u[1]-a[1]; ua[2] = u[2]-a[2];
	vu[0] = v[0]-u[0]; vu[1] = v[1]-u[1]; vu[2] = v[2]-u[2];
	n[0] = ba[1]*vu[2]-ba[2]*vu[1];
	n[1] = ba[2]*vu[0]-ba[0]*vu[2];
	n[2] = ba[0]*vu[1]-ba[1]*vu[0];
	nxua[0] = n[1]*ua[2]-n[2]*ua[1];
	nxua[1] = n[2]*ua[0]-n[0]*ua[2];
	nxua[2] = n[0]*ua[1]-n[1]*ua[0];
	nn = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
	nua = n[0]*ua[0]+n[1]*ua[1]+n[2]*ua[2];

	l1 = vu[0]*nxua[0]+vu[1]*nxua[1]+vu[2]*nxua[2];
	l2 = ba[0]*nxua[0]+ba[1]*nxua[1]+ba[2]*nxua[2];
	if ((l1 > 0.0) & (l1 < nn) & (l2 > 0.0) & (l2 < nn)) return nua*nua/nn;
	return 1e10;
}

//float vtDist2(float p[],float a[],float b[],float c[]) {
float vtDist2( const Vector3<float>& p, const Vector3<float>& a, const Vector3<float>& b, const Vector3<float>& c )
{
	Vector3<float> ba,ca,pa,n,nxpa;
	float nn,npa,l1,l2;

	ba[0] = b[0]-a[0]; ba[1] = b[1]-a[1]; ba[2] = b[2]-a[2];
	ca[0] = c[0]-a[0]; ca[1] = c[1]-a[1]; ca[2] = c[2]-a[2];
	pa[0] = p[0]-a[0]; pa[1] = p[1]-a[1]; pa[2] = p[2]-a[2];
	n[0] = ba[1]*ca[2]-ba[2]*ca[1];
	n[1] = ba[2]*ca[0]-ba[0]*ca[2];
	n[2] = ba[0]*ca[1]-ba[1]*ca[0];
	nn = n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
	npa = n[0]*pa[0]+n[1]*pa[1]+n[2]*pa[2];
	nxpa[0] = n[1]*pa[2]-n[2]*pa[1];
	nxpa[1] = n[2]*pa[0]-n[0]*pa[2];
	nxpa[2] = n[0]*pa[1]-n[1]*pa[0];
	l1 = -ba[0]*nxpa[0]-ba[1]*nxpa[1]-ba[2]*nxpa[2];
	l2 = ca[0]*nxpa[0]+ca[1]*nxpa[1]+ca[2]*nxpa[2];
	if ((l1 > 0.0) & (l1 < nn) & (l2 > 0.0) & (l2 < nn) & (l1+l2 < nn))
		return npa*npa/nn;
	return 1e10;
}



float teDistance2(const leprechaun::geometry::Triangle<float>& T, const Vector3<float>& p, const Vector3<float>& q)
{
	float dmin = 1e10;

	const Vector3<float> a = T[0];
	const Vector3<float> b = T[1];
	const Vector3<float> c = T[2];

	if( teIntersect(T,p,q) ) return 0.0;

	dmin = fminf(vvDist2(a,p),dmin);
	dmin = fminf(vvDist2(a,q),dmin);
	dmin = fminf(vvDist2(b,p),dmin);
	dmin = fminf(vvDist2(b,q),dmin);
	dmin = fminf(vvDist2(c,p),dmin);
	dmin = fminf(vvDist2(c,q),dmin);

	dmin = fminf(veDist2(a,p,q),dmin);
	dmin = fminf(veDist2(b,p,q),dmin);
	dmin = fminf(veDist2(c,p,q),dmin);

	dmin = fminf(veDist2(p,a,b),dmin);
	dmin = fminf(veDist2(p,b,c),dmin);
	dmin = fminf(veDist2(p,c,a),dmin);
	dmin = fminf(veDist2(q,a,b),dmin);
	dmin = fminf(veDist2(q,b,c),dmin);
	dmin = fminf(veDist2(q,c,a),dmin);

	dmin = fminf(eeDist2(a,b,p,q),dmin);
	dmin = fminf(eeDist2(b,c,p,q),dmin);
	dmin = fminf(eeDist2(c,a,p,q),dmin);

	dmin = fminf(vtDist2(p,a,b,c),dmin);
	dmin = fminf(vtDist2(q,a,b,c),dmin);

	return dmin;
}




inline float distPointLine2(const Vector3<float>& p, const Vector3<float>& q,Vector3<float>& c )
{
	return std::min( std::min(vvDist2(p,c),vvDist2(q,c)),veDist2(c,p,q) );
}




} // namespace SV
#endif //SV_GEOMETRY_H
