#ifndef SURFACETRIMMER_INCLUDED
#define SURFACETRIMMER_INCLUDED
#undef ARRAY_DEBUG
#define DIMENSION 3

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "Mesh/PoissonRecon/FEMTree.h"

#include "Mesh/PoissonRecon/MyMiscellany.h"

#include "Mesh/PoissonRecon/CmdLineParser.h"

#include "Mesh/PoissonRecon/MAT.h"

#include "Mesh/PoissonRecon/Geometry.h"

#include "Mesh/PoissonRecon/Ply.h"

#include "Mesh/PoissonRecon/PointStreamData.h"

class SurfaceTrimmer
{
 public:
  MessageWriter messageWriter;

  cmdLineParameter<char*> In{"in"}, Out{"out"};
  cmdLineParameter<int> Smooth{"smooth", 5};
  cmdLineParameter<float> Trim{"trim"}, IslandAreaRatio{"aRatio", 0.001f};
  cmdLineReadable PolygonMesh{"polygonMesh"}, Verbose{"verbose"};

  cmdLineReadable* params[8] = {&In,
                                &Out,
                                &Trim,
                                &PolygonMesh,
                                &Smooth,
                                &IslandAreaRatio,
                                &Verbose,
                                NULL};

  void ShowUsage(char* ex)
  {
    printf("Usage: %s\n", ex);
    printf("\t --%s <input polygon mesh>\n", In.name);
    printf("\t --%s <trimming value>\n", Trim.name);
    printf("\t[--%s <ouput polygon mesh>]\n", Out.name);
    printf("\t[--%s <smoothing iterations>=%d]\n", Smooth.name, Smooth.value);
    printf("\t[--%s <relative area of islands>=%f]\n",
           IslandAreaRatio.name,
           IslandAreaRatio.value);
    printf("\t[--%s]\n", PolygonMesh.name);
    printf("\t[--%s]\n", Verbose.name);
  }

  long long EdgeKey(int key1, int key2)
  {
    if (key1 < key2)
      return (((long long)key1) << 32) | ((long long)key2);
    else
      return (((long long)key2) << 32) | ((long long)key1);
  }

  template <typename Real, typename... VertexData>
  PlyVertexWithData<
      float,
      DIMENSION,
      MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>
  InterpolateVertices(
      const PlyVertexWithData<
          float,
          DIMENSION,
          MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>&
          v1,
      const PlyVertexWithData<
          float,
          DIMENSION,
          MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>&
          v2,
      Real value)
  {
    if (std::get<0>(v1.data.data).data == std::get<0>(v2.data.data).data)
      return (v1 + v2) / Real(2.);
    Real dx = (std::get<0>(v1.data.data).data - value) /
              (std::get<0>(v1.data.data).data - std::get<0>(v2.data.data).data);
    return v1 * (1.f - dx) + v2 * dx;
  }
  template <typename Real, typename... VertexData>
  void SmoothValues(
      std::vector<PlyVertexWithData<
          float,
          DIMENSION,
          MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>>&
          vertices,
      const std::vector<std::vector<int>>& polygons)
  {
    std::vector<int> count(vertices.size());
    std::vector<Real> sums(vertices.size(), 0);
    for (size_t i = 0; i < polygons.size(); i++)
    {
      int sz = int(polygons[i].size());
      for (int j = 0; j < sz; j++)
      {
        int j1 = j, j2 = (j + 1) % sz;
        int v1 = polygons[i][j1], v2 = polygons[i][j2];
        count[v1]++, count[v2]++;
        sums[v1] += std::get<0>(vertices[v2].data.data).data,
            sums[v2] += std::get<0>(vertices[v1].data.data).data;
      }
    }
    for (size_t i = 0; i < vertices.size(); i++)
      std::get<0>(vertices[i].data.data).data =
          (sums[i] + std::get<0>(vertices[i].data.data).data) / (count[i] + 1);
  }
  template <class Real, typename... VertexData>
  void SplitPolygon(
      const std::vector<int>& polygon,
      std::vector<PlyVertexWithData<
          float,
          DIMENSION,
          MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>>&
          vertices,
      std::vector<std::vector<int>>* ltPolygons,
      std::vector<std::vector<int>>* gtPolygons,
      std::vector<bool>* ltFlags,
      std::vector<bool>* gtFlags,
      std::unordered_map<long long, int>& vertexTable,
      Real trimValue)
  {
    int sz = int(polygon.size());
    std::vector<bool> gt(sz);
    int gtCount = 0;
    for (int j = 0; j < sz; j++)
    {
      gt[j] = (std::get<0>(vertices[polygon[j]].data.data).data > trimValue);
      if (gt[j]) gtCount++;
    }
    if (gtCount == sz)
    {
      if (gtPolygons) gtPolygons->push_back(polygon);
      if (gtFlags) gtFlags->push_back(false);
    }
    else if (gtCount == 0)
    {
      if (ltPolygons) ltPolygons->push_back(polygon);
      if (ltFlags) ltFlags->push_back(false);
    }
    else
    {
      int start;
      for (start = 0; start < sz; start++)
        if (gt[start] && !gt[(start + sz - 1) % sz]) break;

      bool gtFlag = true;
      std::vector<int> poly;

      // Add the initial vertex
      {
        int j1 = (start + int(sz) - 1) % sz, j2 = start;
        int v1 = polygon[j1], v2 = polygon[j2];
        int vIdx;
        std::unordered_map<long long, int>::iterator iter =
            vertexTable.find(EdgeKey(v1, v2));
        if (iter == vertexTable.end())
        {
          vertexTable[EdgeKey(v1, v2)] = vIdx = int(vertices.size());
          vertices.push_back(
              InterpolateVertices(vertices[v1], vertices[v2], trimValue));
        }
        else
          vIdx = iter->second;
        poly.push_back(vIdx);
      }

      for (int _j = 0; _j <= sz; _j++)
      {
        int j1 = (_j + start + sz - 1) % sz, j2 = (_j + start) % sz;
        int v1 = polygon[j1], v2 = polygon[j2];
        if (gt[j2] == gtFlag)
          poly.push_back(v2);
        else
        {
          int vIdx;
          std::unordered_map<long long, int>::iterator iter =
              vertexTable.find(EdgeKey(v1, v2));
          if (iter == vertexTable.end())
          {
            vertexTable[EdgeKey(v1, v2)] = vIdx = int(vertices.size());
            vertices.push_back(
                InterpolateVertices(vertices[v1], vertices[v2], trimValue));
          }
          else
            vIdx = iter->second;
          poly.push_back(vIdx);
          if (gtFlag)
          {
            if (gtPolygons) gtPolygons->push_back(poly);
            if (ltFlags) ltFlags->push_back(true);
          }
          else
          {
            if (ltPolygons) ltPolygons->push_back(poly);
            if (gtFlags) gtFlags->push_back(true);
          }
          poly.clear(), poly.push_back(vIdx), poly.push_back(v2);
          gtFlag = !gtFlag;
        }
      }
    }
  }
  template <class Real, class Vertex>
  void Triangulate(const std::vector<Vertex>& vertices,
                   const std::vector<std::vector<int>>& polygons,
                   std::vector<std::vector<int>>& triangles)
  {
    triangles.clear();
    for (size_t i = 0; i < polygons.size(); i++)
      if (polygons.size() > 3)
      {
        std::vector<Point<Real, DIMENSION>> _vertices(polygons[i].size());
        for (int j = 0; j < int(polygons[i].size()); j++)
          _vertices[j] = vertices[polygons[i][j]].point;
        std::vector<TriangleIndex> _triangles =
            MinimalAreaTriangulation<Real, DIMENSION>(
                (ConstPointer(Point<Real, DIMENSION>))GetPointer(_vertices),
                _vertices.size());

        // Add the triangles to the mesh
        size_t idx = triangles.size();
        triangles.resize(idx + _triangles.size());
        for (int j = 0; j < int(_triangles.size()); j++)
        {
          triangles[idx + j].resize(3);
          for (int k = 0; k < 3; k++)
            triangles[idx + j][k] = polygons[i][_triangles[j].idx[k]];
        }
      }
      else if (polygons[i].size() == 3)
        triangles.push_back(polygons[i]);
  }
  template <class Real, class Vertex>
  double PolygonArea(const std::vector<Vertex>& vertices,
                     const std::vector<int>& polygon)
  {
    if (polygon.size() < 3)
      return 0.;
    else if (polygon.size() == 3)
      return Area(vertices[polygon[0]].point,
                  vertices[polygon[1]].point,
                  vertices[polygon[2]].point);
    else
    {
      Point<Real, DIMENSION> center;
      for (size_t i = 0; i < polygon.size(); i++)
        center += vertices[polygon[i]].point;
      center /= Real(polygon.size());
      double area = 0;
      for (size_t i = 0; i < polygon.size(); i++)
        area += Area(center,
                     vertices[polygon[i]].point,
                     vertices[polygon[(i + 1) % polygon.size()]].point);
      return area;
    }
  }

  template <class Vertex>
  void RemoveHangingVertices(std::vector<Vertex>& vertices,
                             std::vector<std::vector<int>>& polygons)
  {
    std::unordered_map<int, int> vMap;
    std::vector<bool> vertexFlags(vertices.size(), false);
    for (size_t i = 0; i < polygons.size(); i++)
      for (size_t j = 0; j < polygons[i].size(); j++)
        vertexFlags[polygons[i][j]] = true;
    int vCount = 0;
    for (int i = 0; i < int(vertices.size()); i++)
      if (vertexFlags[i]) vMap[i] = vCount++;
    for (size_t i = 0; i < polygons.size(); i++)
      for (size_t j = 0; j < polygons[i].size(); j++)
        polygons[i][j] = vMap[polygons[i][j]];

    std::vector<Vertex> _vertices(vCount);
    for (int i = 0; i < int(vertices.size()); i++)
      if (vertexFlags[i]) _vertices[vMap[i]] = vertices[i];
    vertices = _vertices;
  }
  void SetConnectedComponents(const std::vector<std::vector<int>>& polygons,
                              std::vector<std::vector<int>>& components)
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
        std::unordered_map<long long, int>::iterator iter =
            edgeTable.find(eKey);
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
      int p = int(i);
      while (polygonRoots[p] != p) p = polygonRoots[p];
      int root = p;
      p        = int(i);
      while (polygonRoots[p] != p)
      {
        int temp        = polygonRoots[p];
        polygonRoots[p] = root;
        p               = temp;
      }
    }
    int cCount = 0;
    std::unordered_map<int, int> vMap;
    for (int i = 0; i < int(polygonRoots.size()); i++)
      if (polygonRoots[i] == i) vMap[i] = cCount++;
    components.resize(cCount);
    for (int i = 0; i < int(polygonRoots.size()); i++)
      components[vMap[polygonRoots[i]]].push_back(i);
  }

  ///////////////////////////////////////////////////////////
  template <typename... VertexData>
  int Execute(std::vector<double>& vertz_iuv,
              std::vector<double>& clrz_iuv,
              std::vector<int>& tri_indices)
  {
    typedef PlyVertexWithData<
        float,
        DIMENSION,
        MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>
        Vertex;
    float min, max;
    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> polygons;

    int ft, commentNum = 0;
    std::vector<char*> comments;
    char** _comments;
    PlyReadPolygons<Vertex>(In.value,
                            vertices,
                            polygons,
                            Vertex::PlyReadProperties(),
                            Vertex::PlyReadNum,
                            ft,
                            &_comments,
                            &commentNum);
    for (int i = 0; i < commentNum; i++) comments.push_back(_comments[i]);
    for (int i = 0; i < Smooth.value; i++)
      SmoothValues<float>(vertices, polygons);
    min = max = std::get<0>(vertices[0].data.data).data;
    for (size_t i = 0; i < vertices.size(); i++)
      min = std::min<float>(min, std::get<0>(vertices[0].data.data).data),
      max = std::max<float>(max, std::get<0>(vertices[0].data.data).data);
    if (Verbose.set) printf("Value Range: [%f,%f]\n", min, max);

    std::unordered_map<long long, int> vertexTable;
    std::vector<std::vector<int>> ltPolygons, gtPolygons;
    std::vector<bool> ltFlags, gtFlags;

    messageWriter(comments, "*********************************************\n");
    messageWriter(comments, "*********************************************\n");
    messageWriter(
        comments, "** Running Surface Trimmer (Version %s) **\n", VERSION);
    messageWriter(comments, "*********************************************\n");
    messageWriter(comments, "*********************************************\n");
    char str[1024];
    for (int i = 0; params[i]; i++)
      if (params[i]->set)
      {
        params[i]->writeValue(str);
        if (strlen(str))
          messageWriter(comments, "\t--%s %s\n", params[i]->name, str);
        else
          messageWriter(comments, "\t--%s\n", params[i]->name);
      }

    double t = Time();
    for (size_t i = 0; i < polygons.size(); i++)
      SplitPolygon(polygons[i],
                   vertices,
                   &ltPolygons,
                   &gtPolygons,
                   &ltFlags,
                   &gtFlags,
                   vertexTable,
                   Trim.value);
    if (IslandAreaRatio.value > 0)
    {
      std::vector<std::vector<int>> _ltPolygons, _gtPolygons;
      std::vector<std::vector<int>> ltComponents, gtComponents;
      SetConnectedComponents(ltPolygons, ltComponents);
      SetConnectedComponents(gtPolygons, gtComponents);
      std::vector<double> ltAreas(ltComponents.size(), 0.),
          gtAreas(gtComponents.size(), 0.);
      std::vector<bool> ltComponentFlags(ltComponents.size(), false),
          gtComponentFlags(gtComponents.size(), false);
      double area = 0.;
      for (size_t i = 0; i < ltComponents.size(); i++)
      {
        for (size_t j = 0; j < ltComponents[i].size(); j++)
        {
          ltAreas[i] += PolygonArea<float, Vertex>(
              vertices, ltPolygons[ltComponents[i][j]]);
          ltComponentFlags[i] =
              (ltComponentFlags[i] || ltFlags[ltComponents[i][j]]);
        }
        area += ltAreas[i];
      }
      for (size_t i = 0; i < gtComponents.size(); i++)
      {
        for (size_t j = 0; j < gtComponents[i].size(); j++)
        {
          gtAreas[i] += PolygonArea<float, Vertex>(
              vertices, gtPolygons[gtComponents[i][j]]);
          gtComponentFlags[i] =
              (gtComponentFlags[i] || gtFlags[gtComponents[i][j]]);
        }
        area += gtAreas[i];
      }
      for (size_t i = 0; i < ltComponents.size(); i++)
      {
        if (ltAreas[i] < area * IslandAreaRatio.value && ltComponentFlags[i])
          for (size_t j = 0; j < ltComponents[i].size(); j++)
            _gtPolygons.push_back(ltPolygons[ltComponents[i][j]]);
        else
          for (size_t j = 0; j < ltComponents[i].size(); j++)
            _ltPolygons.push_back(ltPolygons[ltComponents[i][j]]);
      }
      for (size_t i = 0; i < gtComponents.size(); i++)
      {
        if (gtAreas[i] < area * IslandAreaRatio.value && gtComponentFlags[i])
          for (size_t j = 0; j < gtComponents[i].size(); j++)
            _ltPolygons.push_back(gtPolygons[gtComponents[i][j]]);
        else
          for (size_t j = 0; j < gtComponents[i].size(); j++)
            _gtPolygons.push_back(gtPolygons[gtComponents[i][j]]);
      }
      ltPolygons = _ltPolygons, gtPolygons = _gtPolygons;
    }
    if (!PolygonMesh.set)
    {
      {
        std::vector<std::vector<int>> polys = ltPolygons;
        Triangulate<float, Vertex>(vertices, ltPolygons, polys),
            ltPolygons = polys;
      }
      {
        std::vector<std::vector<int>> polys = gtPolygons;
        Triangulate<float, Vertex>(vertices, gtPolygons, polys),
            gtPolygons = polys;
      }
    }

    RemoveHangingVertices(vertices, gtPolygons);
    sprintf(comments[commentNum++], "#Trimmed In: %9.1f (s)", Time() - t);

    vertz_iuv.resize(vertices.size() * 3);
    clrz_iuv.resize(vertices.size() * 3);

    for (size_t i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < 3; j++)

      {
        vertz_iuv[i * 3 + j] = vertices[i].point[j];
        clrz_iuv[i * 3 + j]  = std::get<1>(vertices[i].data.data).data[j];
      }
    }

    // std::vector<int>tri_indices;

    for (size_t i = 0; i < gtPolygons.size(); i++)
    {
      size_t size = gtPolygons[i].size();

      for (size_t j = 0; j < size; j++)
      {
        tri_indices.push_back(gtPolygons[i][j]);
      }
    }
    if (Out.set)
      PlyWritePolygons<Vertex>(Out.value,
                               vertices,
                               gtPolygons,
                               Vertex::PlyWriteProperties(),
                               Vertex::PlyWriteNum,
                               ft,
                               &comments[0],
                               (int)comments.size());

    return EXIT_SUCCESS;
  }

  ///////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  template <typename... VertexData>
  int Execute(std::vector<double>& vertz_iuv,
              std::vector<double>& nrmlz_iuv,
              std::vector<double>& clrz_iuv,
              std::vector<double>& roughness_iuv,
              std::vector<double>& metallic_iuv,
              std::vector<int>& tri_indices)
  {
    typedef PlyVertexWithData<
        float,
        DIMENSION,
        MultiPointStreamData<float, PointStreamValue<float>, VertexData...>>
        Vertex;
    float min, max;
    std::vector<Vertex> vertices;
    std::vector<std::vector<int>> polygons;

    int ft, commentNum = 0;
    std::vector<char*> comments;
    char** _comments;
    PlyReadPolygons<Vertex>(In.value,
                            vertices,
                            polygons,
                            Vertex::PlyReadProperties(),
                            Vertex::PlyReadNum,
                            ft,
                            &_comments,
                            &commentNum);
    for (int i = 0; i < commentNum; i++) comments.push_back(_comments[i]);
    for (int i = 0; i < Smooth.value; i++)
      SmoothValues<float>(vertices, polygons);
    min = max = std::get<0>(vertices[0].data.data).data;
    for (size_t i = 0; i < vertices.size(); i++)
      min = std::min<float>(min, std::get<0>(vertices[0].data.data).data),
      max = std::max<float>(max, std::get<0>(vertices[0].data.data).data);
    if (Verbose.set) printf("Value Range: [%f,%f]\n", min, max);

    std::unordered_map<long long, int> vertexTable;
    std::vector<std::vector<int>> ltPolygons, gtPolygons;
    std::vector<bool> ltFlags, gtFlags;

    messageWriter(comments, "*********************************************\n");
    messageWriter(comments, "*********************************************\n");
    messageWriter(
        comments, "** Running Surface Trimmer (Version %s) **\n", VERSION);
    messageWriter(comments, "*********************************************\n");
    messageWriter(comments, "*********************************************\n");
    char str[1024];
    for (int i = 0; params[i]; i++)
      if (params[i]->set)
      {
        params[i]->writeValue(str);
        if (strlen(str))
          messageWriter(comments, "\t--%s %s\n", params[i]->name, str);
        else
          messageWriter(comments, "\t--%s\n", params[i]->name);
      }

    double t = Time();
    for (size_t i = 0; i < polygons.size(); i++)
      SplitPolygon(polygons[i],
                   vertices,
                   &ltPolygons,
                   &gtPolygons,
                   &ltFlags,
                   &gtFlags,
                   vertexTable,
                   Trim.value);
    if (IslandAreaRatio.value > 0)
    {
      std::vector<std::vector<int>> _ltPolygons, _gtPolygons;
      std::vector<std::vector<int>> ltComponents, gtComponents;
      SetConnectedComponents(ltPolygons, ltComponents);
      SetConnectedComponents(gtPolygons, gtComponents);
      std::vector<double> ltAreas(ltComponents.size(), 0.),
          gtAreas(gtComponents.size(), 0.);
      std::vector<bool> ltComponentFlags(ltComponents.size(), false),
          gtComponentFlags(gtComponents.size(), false);
      double area = 0.;
      for (size_t i = 0; i < ltComponents.size(); i++)
      {
        for (size_t j = 0; j < ltComponents[i].size(); j++)
        {
          ltAreas[i] += PolygonArea<float, Vertex>(
              vertices, ltPolygons[ltComponents[i][j]]);
          ltComponentFlags[i] =
              (ltComponentFlags[i] || ltFlags[ltComponents[i][j]]);
        }
        area += ltAreas[i];
      }
      for (size_t i = 0; i < gtComponents.size(); i++)
      {
        for (size_t j = 0; j < gtComponents[i].size(); j++)
        {
          gtAreas[i] += PolygonArea<float, Vertex>(
              vertices, gtPolygons[gtComponents[i][j]]);
          gtComponentFlags[i] =
              (gtComponentFlags[i] || gtFlags[gtComponents[i][j]]);
        }
        area += gtAreas[i];
      }
      for (size_t i = 0; i < ltComponents.size(); i++)
      {
        if (ltAreas[i] < area * IslandAreaRatio.value && ltComponentFlags[i])
          for (size_t j = 0; j < ltComponents[i].size(); j++)
            _gtPolygons.push_back(ltPolygons[ltComponents[i][j]]);
        else
          for (size_t j = 0; j < ltComponents[i].size(); j++)
            _ltPolygons.push_back(ltPolygons[ltComponents[i][j]]);
      }
      for (size_t i = 0; i < gtComponents.size(); i++)
      {
        if (gtAreas[i] < area * IslandAreaRatio.value && gtComponentFlags[i])
          for (size_t j = 0; j < gtComponents[i].size(); j++)
            _ltPolygons.push_back(gtPolygons[gtComponents[i][j]]);
        else
          for (size_t j = 0; j < gtComponents[i].size(); j++)
            _gtPolygons.push_back(gtPolygons[gtComponents[i][j]]);
      }
      ltPolygons = _ltPolygons, gtPolygons = _gtPolygons;
    }
    if (!PolygonMesh.set)
    {
      {
        std::vector<std::vector<int>> polys = ltPolygons;
        Triangulate<float, Vertex>(vertices, ltPolygons, polys),
            ltPolygons = polys;
      }
      {
        std::vector<std::vector<int>> polys = gtPolygons;
        Triangulate<float, Vertex>(vertices, gtPolygons, polys),
            gtPolygons = polys;
      }
    }

    RemoveHangingVertices(vertices, gtPolygons);
    sprintf(comments[commentNum++], "#Trimmed In: %9.1f (s)", Time() - t);

    vertz_iuv.resize(vertices.size() * 3);
    clrz_iuv.resize(vertices.size() * 3);
    nrmlz_iuv.resize(vertices.size() * 3);
    roughness_iuv.resize(vertices.size());
    metallic_iuv.resize(vertices.size());

    for (size_t i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < 3; j++)

      {
        vertz_iuv[i * 3 + j] = vertices[i].point[j];

        nrmlz_iuv[i * 3 + j] = std::get<1>(vertices[i].data.data).data[j];
        clrz_iuv[i * 3 + j]  = std::get<2>(vertices[i].data.data).data[j];
      }
      roughness_iuv[i] = std::get<3>(vertices[i].data.data).data;
      metallic_iuv[i]  = std::get<4>(vertices[i].data.data).data;
    }

    // std::vector<int>tri_indices;

    for (size_t i = 0; i < gtPolygons.size(); i++)
    {
      size_t size = gtPolygons[i].size();

      for (size_t j = 0; j < size; j++)
      {
        tri_indices.push_back(gtPolygons[i][j]);
      }
    }
    if (Out.set)
      PlyWritePolygons<Vertex>(Out.value,
                               vertices,
                               gtPolygons,
                               Vertex::PlyWriteProperties(),
                               Vertex::PlyWriteNum,
                               ft,
                               &comments[0],
                               (int)comments.size());

    return EXIT_SUCCESS;
  }

  int trim(int numOptions,
           char* options[],
           std::vector<double>& vertz_iuv,
           std::vector<double>& clrz_iuv,
           std::vector<int>& tri_indices)
  {
    cmdLineParse(numOptions - 1, &options[1], params);
    messageWriter.echoSTDOUT = Verbose.set;

    if (!In.set || !Trim.set)
    {
      ShowUsage(options[0]);
      return EXIT_FAILURE;
    }
    typedef MultiPointStreamData<float,
                                 PointStreamValue<float>,
                                 PointStreamColor<float>>
        VertexData;

    typedef PlyVertexWithData<float, DIMENSION, VertexData> Vertex;
    bool readFlags[Vertex::PlyReadNum];
    if (!PlyReadHeader(In.value,
                       Vertex::PlyReadProperties(),
                       Vertex::PlyReadNum,
                       readFlags))
      fprintf(stderr, "[ERROR] Failed to read ply header: %s\n", In.value),
          exit(0);

    bool hasValue =
        VertexData::ValidPlyReadProperties<0>(readFlags + DIMENSION);
    bool hasColor =
        VertexData::ValidPlyReadProperties<1>(readFlags + DIMENSION);

    if (!hasValue)
      fprintf(stderr, "[ERROR] Ply file does not contain values\n"), exit(0);

    if (hasColor)
      return Execute<PointStreamColor<float>>(vertz_iuv, clrz_iuv, tri_indices);
  }

  int trim(int numOptions,
           char* options[],
           std::vector<double>& vertz_iuv,
           std::vector<double>& nrmlz_iuv,
           std::vector<double>& clrz_iuv,
           std::vector<double>& roughness_iuv,
           std::vector<double>& metallic_iuv,
           std::vector<int>& tri_indices)
  {
    cmdLineParse(numOptions - 1, &options[1], params);
    messageWriter.echoSTDOUT = Verbose.set;

    if (!In.set || !Trim.set)
    {
      ShowUsage(options[0]);
      return EXIT_FAILURE;
    }
    typedef MultiPointStreamData<float,
                                 PointStreamValue<float>,
                                 PointStreamNormal<float, DIMENSION>,
                                 PointStreamColor<float>,
                                 PointStreamRoughness<float>,
                                 PointStreamMetal<float>>
        VertexData;

    typedef PlyVertexWithData<float, DIMENSION, VertexData> Vertex;
    bool readFlags[Vertex::PlyReadNum];
    if (!PlyReadHeader(In.value,
                       Vertex::PlyReadProperties(),
                       Vertex::PlyReadNum,
                       readFlags))
      fprintf(stderr, "[ERROR] Failed to read ply header: %s\n", In.value),
          exit(0);

    bool hasValue =
        VertexData::ValidPlyReadProperties<0>(readFlags + DIMENSION);
    bool hasNormal =
        VertexData::ValidPlyReadProperties<1>(readFlags + DIMENSION);
    bool hasColor =
        VertexData::ValidPlyReadProperties<2>(readFlags + DIMENSION);
    bool hasRoughness =
        VertexData::ValidPlyReadProperties<3>(readFlags + DIMENSION);
    bool hasMetallic =
        VertexData::ValidPlyReadProperties<4>(readFlags + DIMENSION);

    if (!hasValue)
      fprintf(stderr, "[ERROR] Ply file does not contain values\n"), exit(0);

    if (hasColor)
      if (hasNormal)
      {
        // if (!hasRoughness && !hasMetallic)
        //   return Execute<PointStreamNormal<float, DIMENSION>,
        //                  PointStreamColor<float> >(
        //       vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
        // if (!hasRoughness && hasMetallic)
        //   return Execute<PointStreamNormal<float, DIMENSION>,
        //                  PointStreamColor<float>,
        //                  PointStreamMetal<float> >(
        //       vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
        // if (hasRoughness && !hasMetallic)
        //   return Execute<PointStreamNormal<float, DIMENSION>,
        //                  PointStreamColor<float>,
        //                  PointStreamRoughness<float> >(
        //       vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
        if (hasRoughness && hasMetallic)
          return Execute<PointStreamNormal<float, DIMENSION>,
                         PointStreamColor<float>,
                         PointStreamRoughness<float>,
                         PointStreamMetal<float>>(vertz_iuv,
                                                  nrmlz_iuv,
                                                  clrz_iuv,
                                                  roughness_iuv,
                                                  metallic_iuv,
                                                  tri_indices);
      }
    //   else
    //   {
    //     if (!hasRoughness && !hasMetallic)
    //       return Execute<PointStreamColor<float> >(
    //           vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
    //     if (!hasRoughness && hasMetallic)
    //       return Execute<PointStreamColor<float>, PointStreamMetal<float> >(
    //           vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
    //     if (hasRoughness && !hasMetallic)
    //       return Execute<PointStreamColor<float>, PointStreamRoughness<float>
    //       >(
    //           vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
    //     if (hasRoughness && hasMetallic)
    //       return Execute<PointStreamColor<float>,
    //                      PointStreamRoughness<float>,
    //                      PointStreamMetal<float> >(
    //           vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
    //   }
    // else if (hasNormal)
    //   return Execute<PointStreamNormal<float, DIMENSION> >(
    //       vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
    // else
    //   return Execute<>(vertz_iuv, clrz_iuv, tri_indices, roughness_iuv);
  }
};

#endif  // SURFACETRIMMER_INCLUDED