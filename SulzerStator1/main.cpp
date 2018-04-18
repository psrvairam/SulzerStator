
/*
 * main.cpp
 *
 *  Created on: Jan 27, 2012
 *      Author: renep
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <deque>
#include <algorithm>
using namespace std;
#include <fstream>
#include <iostream>
#include "Param.h"

#include "point.h"
#include "lines.h"
#include "transform.h"
#include "struct.h"
#include "unstruct.h"

#include "OpenGlutVisu.h"

#include "meshTools.h"


struct BoundaryLayer{
	SPLINE internal;
	deque<POINT> internalPts;
	SPLINE external;
	deque<POINT> externalPts;
	UNSTRUCTMESH mesh;
	deque<POINT> meshtopD;
	SPLINE meshtop;
	deque<POINT> meshbotD;
	SPLINE meshbot;
	deque<POINT> meshleadD;
	deque<POINT> meshtrailD;
	deque<POINT> meshtopbD;
	SPLINE meshtopb;
	deque<POINT> meshbotbD;
	SPLINE meshbotb;
	deque<POINT> meshleadbD;
	SPLINE meshleadb;
	deque<POINT> meshtrailbD;
	SPLINE meshtrailb;
};

class CANTILEVERTURBINE : public MESHTOOLS, public VISUAL
{
public:
  CANTILEVERTURBINE(const char *name) : MESHTOOLS(name) {
	  const char* name2 = name;
  }

public:

    vector<string> csv_read_row(istream &in, char delimiter)
    {
        stringstream ss;
        bool inquotes = false;
        vector<std::string> row;//relying on RVO
        while(in.good())
        {
            char c = in.get();
            if (!inquotes && c=='"') //beginquotechar
            {
                inquotes=true;
            }
            else if (inquotes && c=='"') //quotechar
            {
                if ( in.peek() == '"')//2 consecutive quotes resolve to 1
                {
                    ss << (char)in.get();
                }
                else //endquotechar
                {
                    inquotes=false;
                }
            }
            else if (!inquotes && c==delimiter) //end of field
            {
                row.push_back( ss.str() );
                ss.str("");
            }
            else if (!inquotes && (c=='\r' || c=='\n') )
            {
                if(in.peek()=='\n') { in.get(); }
                row.push_back( ss.str() );
                return row;
            }
            else
            {
                ss << c;
            }
        }
    }
    deque<POINT> read_points_csv(string filename){
  	deque<POINT> points;
  	ifstream in(filename.c_str());
  	if (in.fail()) (cout << "File not found" << endl) && 0;
//  	vector<string> row = csv_read_row(in, '\t');
  	while(in.good())
  	{
       vector<string> row = csv_read_row(in, '\t');
       if (row.size() > 0) {
    	   points.push_back(POINT(atof(row[0].c_str()), atof(row[1].c_str()),0));
       }
    }
  	in.close();
  	return points;
    }
	BoundaryLayer BoundaryLayerMesh(deque<POINT> blade, int nElementZ, double fl, double thickBL){
		deque<POINT> internalPts = blade;
		int nElementX = internalPts.size();
		double firstLength = fl;
		double distrMeshRad[nElementZ];
		firstLengthDistr(firstLength, nElementZ, distrMeshRad);
		STRUCTMESH meshBL(nElementX, nElementZ);
		deque<POINT> externalBLpts;
		for (int i=0; i<nElementX; i++)
		{
		  POINT tan;
		  if    ((i==0) || (i==nElementX-1)){
			  tan = internalPts[1] - internalPts[internalPts.size()-2];
		  }
		  else{
			  tan = internalPts[i+1] - internalPts[i-1];
		  }
		  POINT rad(-tan.y, +tan.x);
		  double mag = sqrt(rad.x*rad.x + rad.y*rad.y);
		  rad.x /= mag;
		  rad.y /= mag;

		  for (int j=0; j<nElementZ; j++)
		  {
			meshBL.mesh[i][j] = internalPts[i] + distrMeshRad[j]*thickBL*rad;
			POINT gus(meshBL.mesh[i][j]) ;
			if ((isnan(gus.x)!=0) || (isnan(gus.y)!=0)){
				firstLengthDistr(firstLength, nElementZ, distrMeshRad);
			}
		  }
		  if(i>0)
		  {
			POINT tmp1(meshBL.mesh[i][nElementZ-1]);
			externalBLpts.push_back(tmp1);
		  }
		}
		internalPts.erase(internalPts.end());

		BoundaryLayer bl;
		UNSTRUCTMESH ubl = (meshBL);


		ubl.removeDoublePoints();
		ubl.removeDoubleFaces();
		ubl.UpdateBoundaries();



		for (int i=0; i<ubl.nodes.size(); i++)
			ubl.nodes[i].neig.clear();
		ubl.findNodeNeighbors();


		bl.mesh = ubl;
		bl.internalPts = internalPts;
		bl.externalPts = externalBLpts;
		SPLINE external(externalBLpts);
		SPLINE internal(blade);
		bl.external = external;
		bl.internal = internal;

		return bl;
	}
  void makeMesh()
  {
    // make sure to re-read the file and to delete copies in visu
	    // make sure to re-read the file and to delete copies in visu
    clearParamMap();
    addParamsFromFile("blade.param");
    clearVisu();

    //---------------------- read all the parameters
	string filename = "clean_points.txt";

    double leadtop_par1 = getDoubleParam("LEADTOP_PAR1");
    double leadtop_par2 = getDoubleParam("LEADTOP_PAR2");

    double leadbot_par1 = getDoubleParam("LEADBOT_PAR1");
    double leadbot_par2 = getDoubleParam("LEADBOT_PAR2");

    double trailtop_par1 = getDoubleParam("TRAILTOP_PAR1");
    double trailtop_par2 = getDoubleParam("TRAILTOP_PAR2");

    double trailbot_par1 = getDoubleParam("TRAILBOT_PAR1");
    double trailbot_par2 = getDoubleParam("TRAILBOT_PAR2");

    int nBlades  = getIntParam("NBLADES");
    double Rmid  = getDoubleParam("RADIUS_MIDSPAN");

    int nleadtop = getIntParam("NLEADTOP");
    int nleadbot = getIntParam("NLEADBOT");
    int ntrailtop = getIntParam("NTRAILTOP");
    int ntrailbot = getIntParam("NTRAILBOT");
    int ninlet = getIntParam("NINLET");
    int noutlet = getIntParam("NOUTLET");
    int ntop = getIntParam("NTOP");
    int nbot = getIntParam("NBOT");
    int nelementz = getIntParam("NELEMENTZ");
    double firstlength = getDoubleParam("FIRSTL");
    double thickBL = getDoubleParam("THICKBL");
    string gambitmesh = getStringParam("GAMBITMESH");
    string su2mesh = getStringParam("SU2MESH");
    string fluentmsh = getStringParam("FLUENTMESH");

    bool makemesh = getBoolParam("MAKEMESH");

    //blade
    deque<POINT> blade = read_points_csv("blade_extracted.txt");
    blade.push_back(blade[0]);
    SPLINE tmp2(blade);
    deque<POINT> bladeFinePts = tmp2.discretize(10000, true,false);
    SPLINE bladeS(bladeFinePts);

//    addToDisplay(bladeS);

    int itmid = 3000;
    int itt = 9860;
    int ibmid = 7000;
    int ilt = 5400;
    POINT topmid = bladeFinePts[itmid];
    POINT trailtip = bladeFinePts[itt];
    POINT botmid = bladeFinePts[ibmid];
    POINT leadtip = bladeFinePts[ilt];
//    addToDisplay(topmid	);
//    addToDisplay(trailtip	);
//    addToDisplay(botmid	);
//    addToDisplay(leadtip	);

    deque<POINT> leadtop;
    for (int i=itmid; i<=ilt;i++ ){
    	leadtop.push_back(bladeFinePts[i]);
    }
    SPLINE leadtopS(leadtop);

    deque<POINT> trailbot;
	for (int i=ibmid; i<=itt;i++ ){
		trailbot.push_back(bladeFinePts[i]);
	}
	SPLINE trailbotS(trailbot);
	deque<POINT> leadbot;
	for (int i=ilt; i<=ibmid;i++ ){
		leadbot.push_front(bladeFinePts[i]);
	}
	SPLINE leadbotS(leadbot);
    deque<POINT> trailtop;
	for (int i=itt; i<bladeFinePts.size();i++ ){
		trailtop.push_front(bladeFinePts[i]);
	}
	for (int i=0; i<=itmid;i++ ){
		trailtop.push_front(bladeFinePts[i]);
	}
	SPLINE trailtopS(trailtop);

//	addToDisplay(leadbotS);
//	addToDisplay(trailbotS);
//	addToDisplay(leadtopS);
//	addToDisplay(trailtopS);

	//discretize blade
	deque<POINT> leadbotD = leadbotS.discretize2(nleadbot,1,leadbot_par1,leadbot_par2, true, false);
	deque<POINT> trailbotD = trailbotS.discretize2(ntrailbot,1,trailbot_par1,trailbot_par2, true, false);
	deque<POINT> leadtopD = leadtopS.discretize2(nleadtop,1,leadtop_par1,leadtop_par2, true, false);
	deque<POINT> trailtopD = trailtopS.discretize2(ntrailtop,1,leadbot_par1,leadbot_par2, true, false);
	deque<POINT> bladeD = leadbotD + trailbotD + trailtopD + leadtopD;


	bladeD.push_back(bladeD[0]);

	//create boundary layer
	BoundaryLayer BL = BoundaryLayerMesh(bladeD, nelementz, firstlength, thickBL);
	SPLINE boundarylayerl(BL.externalPts);
	addToDisplay(boundarylayerl);
	boundarylayerl.writeCtrPts("boundaryBL.txt");
	LINE bladeSpline(bladeD);
	bladeSpline.writeCtrPts("blade.txt");


	//boundary
    double transDis = -2*PI*Rmid/double(nBlades);
    POINT tmp(0,transDis,0);

	//periodic
	deque<POINT> wall1 = read_points_csv("boundary_matlab.txt");
	SPLINE wall1S(wall1);
	SPLINE wall2S = wall1S;
	wall2S.translateLine(tmp);

	//outlet
	SPLINE outletS(wall1[0], wall2S.controlPt[0]);
	SPLINE inletS(wall1[wall1.size()-1], wall2S.controlPt[wall2S.controlPt.size()-1]);
	addToDisplay(wall1S);
	addToDisplay(wall2S);
	addToDisplay(inletS);
	addToDisplay(outletS);
	//discretize boundary
	deque<POINT> inletD = inletS.discretize(ninlet,true, false);
	deque<POINT> outletD = outletS.discretize(noutlet,true, false);
	deque<POINT> topD = wall1S.discretize(ntop,true, false);
	deque<POINT> botD = wall2S.discretize(nbot,true, false);
	deque<POINT> boundary = inletD + topD + outletD + botD;
//
	LINE boundaryl(boundary);
	boundaryl.writeCtrPts("boundary.txt");




	addToDisplay(BL.mesh);
	if (makemesh){
		UNSTRUCTMESH outside;
		outside.readFluentMsh2D(gambitmesh.c_str());
		UNSTRUCTMESH mesh2D = BL.mesh+outside ;

		deque<POINT> meshpts = bladeD;
		deque<POINT> periodicboundbot = botD;
		deque<POINT> periodicboundtop = topD;
		deque<POINT> innerpoints_b = inletD;
		deque<POINT> outerpoints_b = outletD;

		for(int i=mesh2D.nfa_i; i<mesh2D.faces.size(); i++)
		{
			if(strcmp (mesh2D.faces[i].name,"noname") == 0)
			{
			  POINT node0= mesh2D.nodes[mesh2D.faces[i].node[0]].pt;
			  POINT node1= mesh2D.nodes[mesh2D.faces[i].node[1]].pt;
			  double tmp_x1= node0.x;
			  double tmp_x2= node1.x;
			  double tmp_y1= node0.y;
			  double tmp_y2= node1.y;

			  if ((fabs(node0.x - inletD[0].x)<1e-8) && (fabs(node1.x - inletD[0].x)<1e-8)) strcpy (mesh2D.faces[i].name,"outlet");
			  else if ((fabs(node0.x - outletD[0].x)<1e-6) && (fabs(node1.x - outletD[0].x)<1e-6)) strcpy (mesh2D.faces[i].name,"inlet");

			  else if (find(meshpts.begin(), meshpts.end(), node1) != meshpts.end() &&
					find(meshpts.begin(), meshpts.end(), node0) != meshpts.end())
					strcpy (mesh2D.faces[i].name,"blade");
			  else {
				  for (int j=0; j < topD.size()-1;j++){
					  if (  (fabs(node1.x - topD[j].x) < 1e-3) && (fabs(node1.y - topD[j].y) < 1e-3)
						 && (fabs(node0.x - topD[j+1].x) < 1e-3) && (fabs(node0.y - topD[j+1].y) < 1e-3)
						 ){
						  strcpy (mesh2D.faces[i].name,"wall1");
					  }
				  }
				  for (int j=0; j < botD.size()-1;j++){
					  if (  (fabs(node0.x - botD[j].x) < 1e-3) && (fabs(node0.y - botD[j].y) < 1e-3)
						 && (fabs(node1.x - botD[j+1].x) < 1e-3) && (fabs(node1.y - botD[j+1].y) < 1e-3)
						 ){
						  strcpy (mesh2D.faces[i].name,"wall2");
					  }
				  }
			  }
			}
		}

		deque<string> BC_names;
		BC_names.push_back("inlet");
		BC_names.push_back("outlet");
		BC_names.push_back("wall1");
		BC_names.push_back("wall2");
		BC_names.push_back("blade");
		mesh2D.writeSU2(su2mesh.c_str(),2, BC_names);
		BC_names.push_back("fluid");
		mesh2D.writeFluentMsh(fluentmsh.c_str(),2, BC_names);

	}

	bCurrTransform = bTransform = false;
    printHelp();
    flagMakeMesh = false;
  }

};




// ------------------------------------------------------------------
//
// Singleton class for opengl call back functions
//
// ------------------------------------------------------------------
class Singleton : public CANTILEVERTURBINE
{
public: // singleton accessors
  static Singleton* instanciate(const char *name)
  {
    if(instance == NULL)
      instance = new Singleton(name);//mesh);
    return instance;
  }

  static Singleton* getInstance()
  {
    if (instance == NULL)  throw(-1);
    else return instance;
  }

  static void reMakeMesh()
  {
    if (flagMakeMesh == true) getInstance()->makeMesh();
    glutPostRedisplay();
  }

  static void display()
  {
    getInstance()->draw();
  }

private:
  static Singleton* instance;  // static instance of EDDY

protected:
  Singleton(const char *name) : CANTILEVERTURBINE(name) {} // overloaded constructor
};

// set singleton instance to zero
Singleton * Singleton::instance = NULL;





int main(int argc, char *argv[])
{
  initOpenGl(argc, argv);

  Singleton *radialComp = Singleton::instanciate(argv[1]);
  radialComp->makeMesh();

  glutDisplayFunc(radialComp->display);
  glutIdleFunc(radialComp->reMakeMesh);
  glutMainLoop();

  return 0;
}






