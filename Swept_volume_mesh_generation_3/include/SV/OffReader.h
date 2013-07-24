// Copyright (c) 2011 Andreas von Dziegielewski and Michael Hemmer (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
//
// Author(s)     : Andreas von Dziegielewski (dziegiel@uni-mainz.de)
//
// ================================================================================

#ifndef LOADOFFFILE_H
#define LOADOFFFILE_H

#include "vecmath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


/// Läd ein OFF-File as_FileName und speichert alle Knoten in ak_Vertices
//* Die Indizierung der Flächen wirdin ak_Indices abgelegt.

void LoadOffFile(const char* as_FileName, std::vector<Vector3d>& ak_Vertices, std::vector<int>& ak_Indices) 
{
  std::cerr << "LoadOffFile(\"" << as_FileName << "\");" << std::endl;

  std::vector<float> lk_Faces;
  int li_TriangleCounter = 0;
  int li_VerticeLength = 0;
  int li_FaceCounter = 0;

	
  std::ifstream lk_InStream(as_FileName);
  if(!lk_InStream)
    {
      std::cerr << "Off-Datei nicht gefunden!" << std::endl;
      return;
    }
  std::string ls_FileHeader;
  std::getline(lk_InStream, ls_FileHeader);
  if (ls_FileHeader.length() > 0 && ls_FileHeader[ls_FileHeader.length()-1] == 0xd) 
    ls_FileHeader.erase (ls_FileHeader.length()-1,1);

  if(ls_FileHeader.compare ("OFF") == 0)
    {
      std::cerr << "Lese OFF-Datei..." << std::endl;
    }
  else 
    {
      std::cerr << "Keine gueltige OFF-Datei!" << std::endl;
      return;
    }

  lk_InStream >> li_VerticeLength;
  std::cerr << "Knoten: " << li_VerticeLength << std::endl;
  lk_InStream >> li_FaceCounter;
  std::cerr << "Flaechen: " << li_FaceCounter << std::endl;
  int li_Edges;
  lk_InStream >> li_Edges;
  std::cerr << "Kanten: " << li_Edges << std::endl;

  for(int i = 0; i < li_VerticeLength; i++) 
    {
      float x, y, z;
      lk_InStream >> x;
      lk_InStream >> y;
      lk_InStream >> z;
      ak_Vertices.push_back(Vector3d(x, y, z));
    }

  int li_Temp1, li_Temp2, li_Temp3;

  for (int i = 0; i < li_FaceCounter; i++) 
    {
      int li_Kanten;
      lk_InStream >> li_Kanten;
      for(int j = 0; j < li_Kanten - 2; j++)
        {
          if (j == 0)
            {
              lk_InStream >> li_Temp1;
              lk_InStream >> li_Temp2;
              lk_InStream >> li_Temp3;
            }
          if (j > 0)
            {
              li_Temp2 = li_Temp3;
              lk_InStream >> li_Temp3;
            }
          li_TriangleCounter++;
          lk_Faces.push_back(float(li_Temp1));
          lk_Faces.push_back(float(li_Temp2));
          lk_Faces.push_back(float(li_Temp3));
        }
    }

  for (int i = 0; i < li_TriangleCounter * 3; i++)
    {
      ak_Indices.push_back(int(lk_Faces[i]));
    }

  //OUT: lk_Vertices, lk_Indices
}

#endif




