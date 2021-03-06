/*
 * point.cpp
 *
 *  Created on: Mar 19, 2012
 *      Author: enrico
 */

#include "point.h"

double POINT::radX() const {return sqrt(y*y + z*z);}
double POINT::radZ() const {return sqrt(x*x + y*y);}
double POINT::phiRadX() const {return atan2(y, z);}
double POINT::phiRadZ() const {return atan2(x, y);}

double POINT::phiDegX() const {return 180.0/PI*atan2(y, z);}
double POINT::phiDegZ() const {return 180.0/PI*atan2(x, y);}

void POINT::rotateDegX(double phi)
{
	phi *= 4.*atan(1.0)/180.0;
	double tmp = y*cos(phi)-z*sin(phi);
	z = y*sin(phi)+z*cos(phi);
	y = tmp;
}

void POINT::rotateDegZ(double phi)
{
	phi *= 4.*atan(1.0)/180.0;
	double tmp = x*cos(phi)-y*sin(phi);
	y = x*sin(phi)+y*cos(phi);
	x = tmp;
}

double POINT::dot(const POINT &p) {return (x*p.x + y*p.y + z*p.z);}

void POINT::drawPoint()
{
	glBegin(GL_POINTS);
	glColor3d(0.2, 0.8, 0.2);
	glVertex3d(x, y, z);
	glEnd();
}

void POINT::print()
{
	printf("xyz: %15.4le\t%15.4le\t%15.4le\n", x, y, z);
}


double mag(const POINT &pt)
{
  return (sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z));
}

POINT crossProd(const POINT &p1, const POINT &p2)
{
  POINT res;
  res.x = p1.y*p2.z - p1.z-p2.y;
  res.y = p1.z*p2.x - p1.x-p2.z;
  res.z = p1.x*p2.y - p1.y-p2.x;
  return res;
}

POINT operator+(const double val, const POINT &pt)
{
  POINT res;
  res.x = val + pt.x;
  res.y = val + pt.y;
  res.z = val + pt.z;
  return res;
}

POINT operator+(const POINT &pt, const double val)
{
  return (val+pt);
}


POINT operator+(const POINT &p1, const POINT &p2)
{
  POINT res;
  res.x = p1.x + p2.x;
  res.y = p1.y + p2.y;
  res.z = p1.z + p2.z;
  return res;
}

POINT operator-(const POINT &p1, const POINT &p2)
{
  POINT res;
  res.x = p1.x - p2.x;
  res.y = p1.y - p2.y;
  res.z = p1.z - p2.z;
  return res;
}

POINT operator-(const double val, const POINT &pt)
{
  POINT res;
  res.x = val - pt.x;
  res.y = val - pt.y;
  res.z = val - pt.z;
  return res;
}

POINT operator-(const POINT &pt, const double val)
{
  return (val-pt);
}



POINT operator+=(POINT &p1, const POINT &p2)
{
  p1.x += p2.x;
  p1.y += p2.y;
  p1.z += p2.z;
  return p1;
}

POINT operator-=(POINT &p1, const POINT &p2)
{
  p1.x -= p2.x;
  p1.y -= p2.y;
  p1.z -= p2.z;
  return p1;
}


bool operator==(POINT &p1, const POINT &p2)
{
  double eps = 1.0e-8;
  return (fabs(p1.x - p2.x) < eps && fabs(p1.y - p2.y) < eps && fabs(p1.z - p2.z) < eps);
}

bool operator!=(POINT &p1, const POINT &p2)
{
  return !(p1 == p2);
}

//bool operator>(const POINT &p1, const POINT &p2)
//{
//  return (p1.x > p2.x);
//}

POINT operator*(const double fact, const POINT &pt)
{
  POINT res;
  res.x = fact*pt.x;
  res.y = fact*pt.y;
  res.z = fact*pt.z;
  return res;
}

POINT operator*(const POINT &pt, const double fact)
{
  return (fact*pt);
}

POINT operator*(const POINT &pt1, const POINT &pt2)
{
  POINT res;
  res.x = pt1.x*pt2.x;
  res.y = pt1.y*pt2.y;
  res.z = pt1.z*pt2.z;
  return res;
}

POINT operator/(const POINT &pt1, const POINT &pt2)
{
  POINT res;

  if (pt2.x != 0.0) res.x = pt1.x/pt2.x;
  if (pt2.y != 0.0) res.y = pt1.y/pt2.y;
  if (pt2.z != 0.0) res.z = pt1.z/pt2.z;
  return res;
}

POINT operator/(const POINT &pt1, const double val)
{
  POINT res;

  if (val != 0.0)
  {
    res.x = pt1.x/val;
    res.y = pt1.y/val;
    res.z = pt1.z/val;
  }
  return res;
}

POINT operator/(const double val, const POINT &pt1)
{
  POINT res;

  if (pt1.x != 0.0) res.x = val/pt1.x;
  if (pt1.y != 0.0) res.y = val/pt1.y;
  if (pt1.z != 0.0) res.z = val/pt1.z;

  return res;
}

POINT rotateDegZ(const double phi, const POINT &p1, const POINT &p2)
{
	POINT p3 = p2-p1;
	p3.rotateDegZ(phi);
  return (p3+p1);
}


deque<POINT> operator+(const deque<POINT> &p1, const deque<POINT> &p2)
{

  deque<POINT> tmp;
  double tmp1 = mag(p1[0]-p2[0]);
  double tmp2 = mag(p1[0]-p2[p2.size()-1]);
  double tmp3 = mag(p1[p1.size()-1]-p2[0]);
  double tmp4= mag(p1[p1.size()-1]-p2[p2.size()-1]);
  double minvalue = min(tmp1,min(tmp2,min(tmp3,tmp4)));

  if (minvalue == tmp1) {
	  for (int i=p1.size()-1; i>=0;i--){
		  tmp.push_back(p1[i]);
	  }
	  if(minvalue < 1e-12) tmp.erase(tmp.end());
	  for (int i=0; i<p2.size();i++){
		  tmp.push_back(p2[i]);
	  }
  }
  else{
	  if(minvalue == tmp2){
		  for (int i=p1.size()-1; i>=0;i--){
			  tmp.push_back(p1[i]);
		  }
		  if(minvalue < 1e-12) tmp.erase(tmp.end());
		  for (int i=p2.size()-1; i>=0;i--){
			  tmp.push_back(p2[i]);
		  }
	  }
	  else {
		  if(minvalue == tmp3){
			  for (int i=0; i<p1.size();i++){
				  tmp.push_back(p1[i]);
			  }
			  if(minvalue < 1e-12 ) tmp.erase(tmp.end());
			  for (int i=0; i<p2.size();i++){
				  tmp.push_back(p2[i]);
			  }
		  }
		  else {
			  for (int i=0; i<p1.size();i++){
				  tmp.push_back(p1[i]);
			  }
			  if(minvalue < 1e-12) tmp.erase(tmp.end());
			  for (int i=p2.size()-1; i>=0;i--){
				  tmp.push_back(p2[i]);
			  }
		  }
	  }
  }

  POINT beginp = tmp[0];
  POINT endp = tmp[tmp.size()-1];
  if (mag(beginp-endp) < 1e-10) tmp.erase(tmp.end());
  return tmp;
}

// read points from file
void readPts(deque<POINT> &pts, const string &name)
{
  FILE *fp = fopen(name.c_str(), "rt");
  if (fp == NULL)
  {
    printf("couldn't open %s\n", name.c_str());
    throw(-1);
  }

  POINT pt;
  while(fscanf(fp, "%lf%lf%lf", &pt.x, &pt.y, &pt.z) != EOF)  { pts.push_back(pt); }
  fclose(fp);
}

