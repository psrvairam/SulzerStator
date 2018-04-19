/*
 * lines.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: enrico
 */
#include "functions.h"
#include "point.h"
#include "lines.h"





//===================================================
//							 LINE CLASS FUNCTIONS
//===================================================
int LINE::size() const { return controlPt.size(); }

void LINE::init() {	calcS(iStart, iEnd); }

POINT & LINE::operator[](int i) { return controlPt[i]; }

POINT const & LINE::operator[](int i) const { return controlPt[i]; }

void LINE::rotateDegX(double phi)
{
	for (int i=0; i<controlPt.size(); i++)
		controlPt[i].rotateDegX(phi);
	//  init();
}

void LINE::rotateDegZ(double phi)
{
	for (int i=0; i<controlPt.size(); i++)
		controlPt[i].rotateDegZ(phi);
}



void LINE::calcS(int i1, int i2)
{
	int N = controlPt.size();

	if (i1 < 0)   i1 = 0;
	if (i2 > N-1) i2 = N-1;

	controlPt[i1].s = 0.0;
	for (int i=i1+1; i<=i2; i++)
		controlPt[i].s = controlPt[i-1].s + mag(controlPt[i]-controlPt[i-1]);
	length = controlPt[i2].s;

	// also assign s for points smaller than i1
	for (int i=i1-1; i>=0; i--)
		controlPt[i].s = controlPt[i+1].s - mag(controlPt[i+1]-controlPt[i]);

	// also assign s for points larger than i2
	for (int i=i2+1; i<N; i++)
		controlPt[i].s = controlPt[i-1].s + mag(controlPt[i]-controlPt[i-1]);

	// normalize
	for (int i=0; i<N; i++)
		controlPt[i].s /= length;
}

POINT LINE::calcPoint(double t, double offset1, double offset2) const
{
	int i=0;

	t = t*(offset2-offset1) + offset1;

	while ((t>=controlPt[i].s) && (i<controlPt.size()-1)) {i++;}

	double fact = (t-controlPt[i-1].s)/(controlPt[i].s-controlPt[i-1].s);
	POINT pt = controlPt[i-1] + fact*(controlPt[i]-controlPt[i-1]);
	return pt;
}


void LINE::split(deque<POINT> &line1, deque<POINT> &line2, double t)
{
	POINT splPt = calcPoint(t);

	int I=1;
	while (controlPt[I].s<t) I++;

	for (int i=0; i<I; i++)
		line1.push_back(controlPt[i]);
	line1.push_back(splPt);

	//	line2.push_back(splPt);     //GUSTAVOOOOOOO 23/06/2016      Line 2 had twice the same points!!!
	for (int i=I; i<controlPt.size(); i++)
		line2.push_back(controlPt[i]);
}

void LINE::splitStephan(LINE intline, deque<POINT> &line1,deque<POINT> &line2)
{
	double t;
	int I;
	POINT intPt = intersectXY(intline, true);
	for (int i=0; i<controlPt.size(); i++){
		if (controlPt[i]==intPt) {I=i; break;}
	}

	for (int i=0; i<=I; i++)
		line1.push_back(controlPt[i]);
	for (int i=I; i<controlPt.size(); i++)
		line2.push_back(controlPt[i]);

	//remove any double points
	unique(line2.begin(), line2.begin());
	deque<POINT>::iterator last = unique(line1.begin(), line1.end());
	line1.resize( distance(line1.begin(),last));
	last = unique(line2.begin(), line2.end());
	line2.resize( distance(line2.begin(),last));
	reverse(line2.begin(), line2.end());
}

void LINE::splitStephan2(LINE intline1,LINE intline2, deque<POINT> &line1,deque<POINT> &line2, int &istart, int &iend)
{
	double s1, s2;
	SPLINE tmps(controlPt);

	POINT intPt1 = tmps.intersectXY(intline1, true);
	POINT intPt2 = tmps.intersectXY(intline2, true);
	for (int i=0; i<tmps.controlPt.size(); i++){
		if (mag(tmps.controlPt[i]-intPt1)<1e-8) s1 = tmps.controlPt[i].s;
		else if (mag(tmps.controlPt[i]-intPt2)<1e-8) s2 = tmps.controlPt[i].s;
	}

	double ss, bs;
	if (s1 < s2) {ss = s1; bs = s2;}
	else {ss=s2; bs=s1;}


	POINT ssplPt1 = calcPoint(ss);
	POINT bsplPt2 = calcPoint(bs);

	int i1, i2;
	for (int i=0; i<controlPt.size(); i++){
		if (controlPt[i].s < ss) i1=i;
		else if (controlPt[i].s > bs) {i2=i; break;}
	}
	cout << i1 << '\t' << i2 << '\n';
	deque<POINT> line1t, line2t;
	for (int i=0; i<=i1;i++)
		line1t.push_back(controlPt[i]);
	for (int i=controlPt.size()-1; i>=i2;i--){
		line1t.push_front(controlPt[i]);
	}
	for (int i=i1+1; i<=i2-1; i++){
		line2t.push_back(controlPt[i]);
	}
	//
	if (line2t.size() < line1t.size()) {line2 = line2t; line1 = line1t;}
	else {line1 = line2t; line2 = line1t;};
	//remove any double points
	unique(line2.begin(), line2.begin());
	deque<POINT>::iterator last = unique(line1.begin(), line1.end());
	line1.resize( distance(line1.begin(),last));
	last = unique(line2.begin(), line2.end());
	line2.resize( distance(line2.begin(),last));
	istart = i1;
	iend = i2;


}
void LINE::splitStephan3(POINT p1, POINT p2, deque<POINT> &line1,deque<POINT> &line2)
{
	double s1, s2;
	int i1, i2;

	POINT intPt1 = p1;
	POINT intPt2 = p2;
	for (int i=0; i<controlPt.size(); i++){
		if (controlPt[i]==intPt1) {
			s1 = controlPt[i].s;
			i1 = i;
		}
		else if (controlPt[i]==intPt2) {
			s2 = controlPt[i].s;
			i2 = i;
		}
	}

	int is, ie;
	if (i1 < i2) {
		is = i1;
		ie = i2;
	}
	else {
		is = i2;
		ie = i1;
	}
	cout << is << '\t' << ie << '\n';


	deque<POINT> line1t, line2t;
	for (int i=0; i<=is;i++)
		line1t.push_back(controlPt[i]);
	for (int i=controlPt.size()-1; i>=ie;i--){
		line1t.push_front(controlPt[i]);
	}
	for (int i=is; i<=ie; i++){
		line2t.push_back(controlPt[i]);
	}
	if (line2t.size() < line1t.size()) {line2 = line2t; line1 = line1t;}
	else {line1 = line2t; line2 = line1t;};
//	//remove any double points
	unique(line2.begin(), line2.begin());
	deque<POINT>::iterator last = unique(line1.begin(), line1.end());
	line1.resize( distance(line1.begin(),last));
	last = unique(line2.begin(), line2.end());
	line2.resize( distance(line2.begin(),last));
}

double LINE::intersect(const POINT &p11, const POINT &p12, const POINT &p21, const POINT &p22)
{
	double denom  = ((p22.y - p21.y)*(p12.x - p11.x)) - ((p22.x - p21.x)*(p12.y - p11.y));
	double nume_a = ((p22.x - p21.x)*(p11.y - p21.y)) - ((p22.y - p21.y)*(p11.x - p21.x));
	double nume_b = ((p12.x - p11.x)*(p11.y - p21.y)) - ((p12.y - p11.y)*(p11.x - p21.x));

	if (denom != 0.0)
	{
		double ua = nume_a/denom;
		double ub = nume_b/denom;
		if (ua >= 0.0 && ua <= 1.0 && ub >= 0.0 && ub <= 1.0)  return ua;
	}
	return -1.0;
}

/*! \brief Computes intersection point of two lines using their control points.
 *
 *  Intersection point is only accurate if both lines are of type LINE
 */
POINT LINE::intersectXY(LINE& line2, bool insertPoint)
{
	POINT intersectPt(0.0, 0.0, 0.0);

	int i = 0, j = 0;
	double ua;
	bool found = false;

	while ((j<line2.size()-1) && !found)
	{
		if ((ua = intersect(controlPt[i], controlPt[i+1], line2.controlPt[j], line2.controlPt[j+1])) != -1.0)   found = true;
		else
		{
			i++;
			if (i>=this->size()-1)
			{
				i = 0;
				j++;
			}
		}
	}

	intersectPt = controlPt[i] + ua*(controlPt[i+1]-controlPt[i]);

	if (insertPoint == true)
	{
		controlPt.insert(controlPt.begin()+i+1, intersectPt);
		iEnd++;
		init();
		//		cout<<"Inserting the point!! "<<endl;
	}

	return intersectPt;
}

/*! \brief Computes intersection point of two lines using their non-dimensional coordinate s.
 *
 *  To be used for SPLINE and BEZIER curves
 *  \param line2 is the line to perform intersection with
 *  \param insertPoint boolean to set if intersection point should be added
 *  \param eps is the truncation error of the iteration to compute the intersection point
 */
POINT LINE::intersectXY(LINE& line2, bool insertPoint, double eps)
{
	POINT intersectPt(0.0, 0.0, 0.0);

	double s1  = 0.0, s1e = 1.0;
	double s1Start = s1;
	double s2  = 0.0, s2e = 1.0;
	double nInterval = 10.0;
	double ds = 1./nInterval;
	double sIntersect, ua;

	while (ds > eps)
	{
		if ((ua=intersect(calcPoint(s1), calcPoint(s1+ds), line2.calcPoint(s2), line2.calcPoint(s2+ds))) != -1.0)
		{
			sIntersect = s1+ua*ds;
			s1e = min(s1+2.0*ds, 1.0);
			s2e = min(s2+2.0*ds, 1.0);
			s1 = s1Start = max(s1-1.*ds, 0.0);
			s2 = max(s2-1.*ds, 0.0);
			ds /= nInterval;
		}
		else
		{
			s1 += ds;
			if (s1>=s1e)
			{
				s1 = s1Start;
				s2 += ds;
				if (s2 >= s2e) return POINT(1.0e99);
			}
		}
	}

	intersectPt = calcPoint(sIntersect);

	if (insertPoint == true)
	{
		int jInsert = 0;
		while(s1 > controlPt[jInsert].s) jInsert++;

		controlPt.insert(controlPt.begin()+jInsert,intersectPt);
		iEnd++;
		init();
	}

	return intersectPt;
}

/*! \brief Computes the orthogonal projection of a point on the line
 *
 *  \param pt is the point to be projected
 */
POINT LINE::project2D(POINT &pt)
{
	double dist = 1e20;
	int ind = 0, indmin = 0;
	deque<POINT> proPts;

	for (int i=1; i<controlPt.size(); i++)
	{
		int inside = 0;

		POINT p1 = controlPt[i-1], p2 = controlPt[i];

		double m1 = (p2.y-p1.y)/(p2.x-p1.x + 1e-10) + 1e-10, m2 = -1/m1;
		double xPro = (p1.y-pt.y+m2*pt.x-m1*p1.x) / (m2-m1);
		double yPro = p1.y + m1*(xPro-p1.x);
		POINT pPro = POINT(xPro,yPro,0.0);

		if (mag(pPro-p1)/mag(p2-p1)<=1)
			inside = 1;

		if (inside==1)
		{
			proPts.push_back(pPro);
			double tmpDist = mag(pPro-pt);

			if (tmpDist<dist)
			{
				dist = tmpDist;
				indmin = ind;
				ind++;
			}
		}
	}

	return (proPts[indmin]);
}
void LINE::translateLine(POINT &pt)
{
	for (int i=0; i<controlPt.size(); i++)
	{
		controlPt[i]+=pt;
	}
}

void LINE::write(const char *name, int N)
{
	FILE *fp = fopen(name, "wt");
	for (int i=0; i<N; i++)
	{
		POINT pt = calcPoint((double)i/(N-1.0));
		fprintf(fp, "%.6lf\t%.6lf\t%.6lf\n", pt.x, pt.y, pt.z);
	}
	fclose(fp);
}

void LINE::write(const char *name, LINE &distr, int N)
{
	FILE *fp = fopen(name, "wt");
	for (int i=0; i<N; i++)
	{
		double t = (double)i/(N-1.0);
		t = distr.calcPoint(t).y;
		POINT pt = calcPoint(t);
		fprintf(fp, "%.6lf\t%.6lf\t%.6lf\n", pt.x, pt.y, pt.z);
	}
	fclose(fp);
}


void LINE::writeCtrPts(const char *name)
{
	FILE *fp = fopen(name, "wt");
	for (int i=0; i<controlPt.size(); i++)
		fprintf(fp, "%.10le\t%.10le\t%.10le\n", controlPt[i].x, controlPt[i].y, controlPt[i].z);
	fclose(fp);
}

void LINE::writeCtrPtsAnsys(const char *name)
{
	FILE *fp = fopen(name, "wt");
//	fprintf(fp, "Group 1\n");
	for (int i=0; i<controlPt.size(); i++)
		fprintf(fp, "%i %i %.10le %.10le %.10le\n",1,i+1, controlPt[i].x, controlPt[i].y, controlPt[i].z);
	fclose(fp);
}


deque<POINT> LINE::discretize(int npoints, bool startend, bool reverse) {
	deque<POINT> points;
	int k, l;
	double s;
	if (startend) { k=0; l=npoints;}
	else { k=1; l=npoints-1;}
	for(int i=k; i<= l; i++) {
		s = (double)i/(double)npoints;
		if(reverse){ points.push_back(calcPoint(s)); }
		else{ points.push_front(calcPoint(s)); }
	}
	return(points);
}
deque<POINT> LINE::discretize2(const int imax,
		const int interp, const double par1, const double par2,  bool startend, bool reversebool) {
	double t[imax];

	for (int i=0; i<imax; i++)
		t[i] = (double)i/(double)(imax-1);


	if (interp==1) // tanh
			{
		for (int i=0; i<imax; i++)
			t[i] = stretchingTanh(t[i], par1, par2);
			}
	else if (interp==2) // atanh
	{
		for (int i=0; i<imax; i++)
			t[i] = stretchingAtanh(t[i], par2, par1);
	}
	else if (interp==3) // first-last length
	{
		firstLastLengthDistr(par1/length, par2/length, imax, t);
	}


	// thickness of the mesh
	LINE line1(imax); // this is the line along the blade geometry
	for (int i=0; i<line1.size(); i++)
	{
		POINT pt = calcPoint(t[i]);
		line1[i] = pt;
	}
	if (!startend) {
		line1.controlPt.erase(line1.controlPt.begin());
		line1.controlPt.erase(line1.controlPt.end());
	}
	if (reversebool) {
		reverse(line1.controlPt.begin(), line1.controlPt.end());
	}

	return(line1.controlPt);
}


void LINE::drawLine()
{
	glLineWidth(3.0);
	glBegin(GL_LINE_STRIP);
	glColor3d(1.0, 0.0, 0.0);
	for (int i=0; i<controlPt.size(); i++)
		glVertex3d(controlPt[i].x, controlPt[i].y, controlPt[i].z);
	glEnd();
	glLineWidth(1.0);
}


//===================================================
//							 BEZIER CLASS FUNCTIONS
//===================================================
POINT BEZIER::calcPoint(double t, double offset1, double offset2) const
{
	t = t*(offset2-offset1) + offset1;
	POINT pt;
	int N = (int)controlPt.size()-1;
	for (int i=0; i<=N; i++)
		pt += (double)binomialCoeff(N, i)*pow(1.-t, double(N-i))*pow(t, double(i))*controlPt[i];
	return pt;
}

void BEZIER::drawLine() // the display of the bezier line is wrong, just connects the points like a line
{
	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	glColor3d(0.0, 1.0, 1.0);
	for (int i=0; i<200; i++)
	{
		double fact = (double)i/199.0;
		glVertex3d(calcPoint(fact).x, calcPoint(fact).y, calcPoint(fact).z);
	}
	glEnd();
	glLineWidth(1.0);
}



//===================================================
//							 SPLINE CLASS FUNCTIONS
//===================================================
void SPLINE::init()
{
	calcS(iStart, iEnd);
	calcDerivative();
}

POINT SPLINE::calcPoint(double t, double offset1, double offset2) const
{
	t = t*(offset2-offset1) + offset1;

	int klo = 0;
	int khi = controlPt.size()-1;

	while (khi-klo > 1)
	{
		int k = (khi+klo) >> 1;
		if (controlPt[k].s > t)  khi = k;
		else                     klo = k;
	}

	double h = controlPt[khi].s - controlPt[klo].s;
	if (h == 0.0) printf("Bad xa input to routine splint\n");

	double a = (controlPt[khi].s - t)/h;
	double b = (t - controlPt[klo].s)/h;
	POINT pt = a*controlPt[klo] + b*controlPt[khi] + (h*h/6.0)*((a*a*a-a)*deriv[klo]+(b*b*b-b)*deriv[khi]);
	return pt;


//    //=====================================================
//    // Extract the derivative of the splined function.
//    // If x0 is outside range return the closest value.
//    //=====================================================
//
//    T dfdx(T x0){
//      int klo, khi;
//      T h, b, a, result;
//
//      // Calculate x[klo] < x0 < x[khi]
//      klo = lower_index(x0);
//      khi = std::min(klo+1,n-1);
//      h = x[khi]-x[klo];
//
//      // Check for error
//      if (h == 0.0){
//        std::cout << "Error in Spline<" << name << "> dfdx(x); h = 0; x-values must be distict!" << std::endl;
//        exit(1);
//      }
//
//      // Calculate interpolation value
//      a = (x[khi]-x0)/h;
//      b = (x0-x[klo])/h;
//      result = (y[khi]-y[klo])/h + h/6.0*(-(3*a*a-1)*y2[klo] + (3*b*b-1)*y2[khi]);
//      return result;
}

void SPLINE::calcDerivative()
{
	int n = controlPt.size();
	deriv.resize(n);
	POINT u[n-1];

	deriv[0] = u[0] = POINT(0.0, 0.0, 0.0);
	for(int i = 1; i < n-1; i++)
	{
		double sig = (controlPt[i].s - controlPt[i-1].s)/(controlPt[i+1].s - controlPt[i-1].s);
		POINT p = sig*deriv[i-1] + 2.0;
		deriv[i] = (sig - 1.0)/p;
		u[i] =  (controlPt[i+1] - controlPt[i])/(controlPt[i+1].s - controlPt[i].s)
        		   -(controlPt[i] - controlPt[i-1])/(controlPt[i].s - controlPt[i-1].s);
		u[i] = (6.0*u[i]/(controlPt[i+1].s - controlPt[i-1].s) - sig*u[i-1])/p;
	}

	POINT qn, un;
	qn = un = POINT(0.0, 0.0, 0.0);
	deriv[n-1] = (un - qn*u[n-2])/(qn*deriv[n-2] + 1.0);
	for(int k = n-2; k >= 0; k--)
		deriv[k] = deriv[k]*deriv[k+1] + u[k];
}

double SPLINE::findS(POINT p)
{
	//====================================Newton Raphson
	double x0, xn;
	double dx = 0.00001;
	double fx, dfx;
	x0=0.0;
	int n=0;
	double tol=100;
	while(n<10000 && tol>0.000000001)
	{
		POINT pts_x0   =  calcPoint(x0);
		POINT pts_x0dx =  calcPoint(x0+dx);
		fx= mag(pts_x0-p);
		dfx= (mag(pts_x0dx-p)-mag(pts_x0-p))/dx;
		xn= x0 - (fx/dfx);
		n++;
		tol= fabs(xn - x0);
		x0=xn;
	}
	//====================================
	if(n<10000 && tol<0.000000001)  return xn;
	else                            {printf("Point is not part of the spline! The tolerance was not achieved (function: SPLINE::findSError) %.10le\n", tol);     return xn;}
}

/*! \brief spline offset
 *
 */
SPLINE SPLINE::offsetRadial(const double thickness)
{
	deque<POINT> newPts;

	for (int i=0; i<controlPt.size(); i++)
		newPts.push_back(controlPt[i]+thickness*calcNorm2D(controlPt[i].s));

	SPLINE newSpl(newPts, iStart, iEnd);

	return (newSpl);

}
POINT SPLINE::calcNorm2D(double t) const
{
	double ds = 1.0e-8; //0.005 ; //0.0001;
	POINT pt1 = calcPoint(t-ds);
	POINT pt2 = calcPoint(t+ds);
	POINT tang = (pt2-pt1)/(2.0*ds);
	double tmp = tang.x;
	tang.x = tang.y;
	tang.y = -tmp;
	return tang/mag(tang);
}
POINT SPLINE::calcTang2D(double t) const
{
	double ds = 1.0e-8; //0.005 ; //0.0001;
	POINT pt1 = calcPoint(t-ds);
	POINT pt2 = calcPoint(t+ds);
	POINT tang = (pt2-pt1)/(2.0*ds);
	return tang/mag(tang);
}


POINT SPLINE::calcNorm2D_new(double t) const
{
	double ds = 0.000001; //0.005 ; //0.0001;
	double angle = atan(calcPoint(t).y/calcPoint(t).x);
	double angle2 = (angle * 180.0)/M_PI;
	POINT pt1 = calcPoint(t-ds);
	POINT pt2 = calcPoint(t+ds);
	POINT tang = (pt2-pt1)/(2.0*ds);
	//POINT tang = news;
	//  tang.rotateDegZ(90);
	double tmp = tang.x;

	if (tang.y >= 0 & tang.x>=0){
		tang.x = tang.y;
		tang.y = tmp;
	}
	else{
		if (tang.y <= 0 & tang.x>=0){
			tang.x = -tang.y;
			tang.y = tmp;
		}
		else {
			if (tang.y <= 0 & tang.x>=0){
				tang.x = -tang.y;
				tang.y = -tmp;
			}
			else{
				tang.x = -tang.y;
				tang.y = tmp;

			}
		}
	}

	return tang/mag(tang);
}



deque<POINT> SPLINE::discretizeSplineS(int npoints, bool startend, bool reverse, double sstart, double send) {
	deque<POINT> points;
	int k, l;
	double s;
	if (startend) {
		k=0;
		l=1;
	} else {
		k=1;
		l=0;
	}
	for(int i=k; i< (npoints+l); i++){
		s = (send-sstart)*(double)i/(double)npoints+sstart;
		if(reverse){
			points.push_back(calcPoint(s));
		}
		else{
			points.push_front(calcPoint(s));
		}
	}
	return(points);
}

void SPLINE::drawLine()
{
	glLineWidth(3.0);
	glBegin(GL_LINE_STRIP);
	glColor3d(1.0, 0.0, 0.5);
	for (int i=0; i<500; i++)
	{
		double fact = (double)i/499.0;
		glVertex3d(calcPoint(fact).x, calcPoint(fact).y, calcPoint(fact).z);
	}
	glEnd();
	glLineWidth(1.0);
}



//===================================================
//							 POLYGON CLASS FUNCTIONS
//===================================================
bool POLYGON::pointInsidePolygon(POINT &pt)
{
	// check if the point pt is inside the polygon
	int intCount=0;

	for (int i=0; i<controlPt.size()-1; i++)
	{
		POINT p1 = controlPt[i], p2 = controlPt[i+1];

		double intX;

		double m = (p2.y-p1.y)/(p2.x-p1.x+1e-20);
		if (m==0) m+=1e-20;

		intX = (pt.y - p1.y)/m + p1.x;

		if ((intX>pt.x) && (intX>=min(p1.x,p2.x)) && (intX<=max(p1.x,p2.x)))
			intCount++;
	}

	bool inside;

	if (intCount%2 == 0)
		inside = false;
	else inside = true;

	return(inside);
}

