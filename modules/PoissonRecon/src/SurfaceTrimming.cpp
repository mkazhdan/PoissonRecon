/*
Copyright (c) 2013, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of
conditions and the following disclaimer. Redistributions in binary form must
reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its
contributors
may be used to endorse or promote products derived from this software without
specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED
WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH
DAMAGE.
*/

#include "Mesh/PoissonRecon/SurfaceTrimming.h"
#include <iostream>

void SurfaceTrimming::ShowUsage(char* ex)
{
  printf("Usage: %s\n", ex);
  printf("\t --%s <input polygon mesh>\n", In.name);
  printf("\t[--%s <trimming value>]\n", Trim.name);
  printf("\t[--%s <ouput polygon mesh>]\n", Out.name);
  printf("\t[--%s <smoothing iterations>=%d]\n", Smooth.name, Smooth.value);
  printf("\t[--%s <relative area of islands>=%f]\n",
         IslandAreaRatio.name,
         IslandAreaRatio.value);
  printf("\t[--%s]\n", PolygonMesh.name);
}

long long SurfaceTrimming::EdgeKey(int key1, int key2)
{
  if (key1 < key2)
    return (((long long)key1) << 32) | ((long long)key2);
  else
    return (((long long)key2) << 32) | ((long long)key1);
}

void SurfaceTrimming::SetConnectedComponents(
    const std::vector<std::vector<int> >& polygons,
    std::vector<std::vector<int> >& components)
{
  std::vector<int> polygonRoots(polygons.size());
  for (size_t i = 0; i < polygons.size(); i++) polygonRoots[i] = int(i);
  std::unordered_map<long long, int> edgeTable;
  for (size_t i = 0; i < polygons.size(); i++)
  {
    int sz = int(polygons[i].size());
    for (int j = 0; j < sz; j++)
    {
      int j1 = j, j2 = (j + 1) % sz;
      int v1 = polygons[i][j1], v2 = polygons[i][j2];
      long long eKey = EdgeKey(v1, v2);
      std::unordered_map<long long, int>::iterator iter = edgeTable.find(eKey);
      if (iter == edgeTable.end())
        edgeTable[eKey] = int(i);
      else
      {
        int p = iter->second;
        while (polygonRoots[p] != p)
        {
          int temp        = polygonRoots[p];
          polygonRoots[p] = int(i);
          p               = temp;
        }
        polygonRoots[p] = int(i);
      }
    }
  }
  for (size_t i = 0; i < polygonRoots.size(); i++)
  {
    int p                          = int(i);
    while (polygonRoots[p] != p) p = polygonRoots[p];
    int root                       = p;
    p                              = int(i);
    while (polygonRoots[p] != p)
    {
      int temp        = polygonRoots[p];
      polygonRoots[p] = root;
      p               = temp;
    }
  }
  int cCount = 0;
  std::unordered_map<int, int> vMap;
  for (int i                          = 0; i < int(polygonRoots.size()); i++)
    if (polygonRoots[i] == i) vMap[i] = cCount++;
  components.resize(cCount);
  for (int i = 0; i < int(polygonRoots.size()); i++)
    components[vMap[polygonRoots[i]]].push_back(i);
}

int SurfaceTrimming::trim(int num_options,
                          char* options[],
                          std::vector<double>& vertz_iuv,
                          std::vector<double>& clrz_iuv,
                          std::vector<double>& tri_indices)
{
  int paramNum = sizeof(params) / sizeof(cmdLineReadable*);
  cmdLineParse(num_options, &options[0], paramNum, params, 0);

  // if (!In.set || !Trim.set)
  if (!Trim.set)
  {
    ShowUsage(options[0]);
    return EXIT_FAILURE;
  }
  bool readFlags[PlyColorAndValueVertex<float>::ReadComponents];
  if (!PlyReadHeader(In.value,
                     PlyColorAndValueVertex<float>::ReadProperties,
                     PlyColorAndValueVertex<float>::ReadComponents,
                     readFlags))
    fprintf(stderr, "[ERROR] Failed to read ply header: %s\n", In.value),
        exit(0);

  bool hasValue = readFlags[3];
  bool hasColor = (readFlags[4] || readFlags[7]) &&
                  (readFlags[5] || readFlags[8]) &&
                  (readFlags[6] || readFlags[9]);

  // std::cout << ";sdkfjh" << std::endl;
  if (!hasValue)
    fprintf(stderr, "[ERROR] Ply file does not contain values\n"), exit(0);
  if (hasColor)
    return Execute<PlyColorAndValueVertex<float> >(
        vertz_iuv, clrz_iuv, tri_indices);
  else
    return Execute<PlyValueVertex<float> >(vertz_iuv, tri_indices);
}
