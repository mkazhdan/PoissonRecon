/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
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
#ifndef POISSONRECON_INCLUDED
#define POISSONRECON_INCLUDED

#undef FAST_COMPILE
#undef ARRAY_DEBUG
#define BRUNO_LEVY_FIX
#define FOR_RELEASE

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#if defined(_WIN32) || defined(_WIN64)
#include <Psapi.h>
#include <Windows.h>
#endif  // _WIN32 || _WIN64
#include "Mesh/PoissonRecon/CmdLineParser.h"
#include "Mesh/PoissonRecon/MarchingCubes.h"
#include "Mesh/PoissonRecon/MemoryUsage.h"
#include "Mesh/PoissonRecon/MyTime.h"
#include "Mesh/PoissonRecon/Octree.h"
#include "Mesh/PoissonRecon/PPolynomial.h"
#include "Mesh/PoissonRecon/Ply.h"
//#include "Mesh/PoissonRecon/PoissonRecon.h"
#include "Mesh/PoissonRecon/SparseMatrix.h"
#ifdef _OPENMP
#include "omp.h"
#endif  // _OPENMP
void DumpOutput(const char* format, ...);
void DumpOutput2(std::vector<char*>& comments, const char* format, ...);
#include "Mesh/PoissonRecon/MultiGridOctreeData.h"

#define DEFAULT_FULL_DEPTH 5

#define XSTR(x) STR(x)
#define STR(x) #x
#if DEFAULT_FULL_DEPTH
#pragma message("[WARNING] Setting default full depth to " XSTR( \
    DEFAULT_FULL_DEPTH))
#endif  // DEFAULT_FULL_DEPTH

#include <stdarg.h>

class PoissonRecon
{
 public:
  // PoissonRecon(){};
  //~PoissonRecon(){};
  char* outputFile = NULL;
  int echoStdout   = 0;
  void DumpOutput(const char* format, ...);

  void DumpOutput2(std::vector<char*>& comments, const char* format, ...);

  cmdLineString In{"in"}, Out{"out"}, TempDir{"tempDir"}, VoxelGrid{"voxel"},
      XForm{"xForm"};

  cmdLineReadable
#if defined(_WIN32) || defined(_WIN64)
      Performance{"performance"},
#endif  // _WIN32 || _WIN64
      ShowResidual{"showResidual"}, NoComments{"noComments"},
      PolygonMesh{"polygonMesh"}, Confidence{"confidence"},
      NormalWeights{"nWeights"}, NonManifold{"nonManifold"}, ASCII{"ascii"},
      Density{"density"}, LinearFit{"linearFit"}, PrimalVoxel{"primalVoxel"},
#ifndef FAST_COMPILE
      Double{"double"},
#endif  // !FAST_COMPILE
      Verbose{"verbose"};

  cmdLineInt
#ifndef FAST_COMPILE
      Degree{"degree", 2},
#endif  // !FAST_COMPILE
      Depth{"depth", 8}, CGDepth{"cgDepth", 0}, KernelDepth{"kernelDepth"},
      AdaptiveExponent{"adaptiveExp", 1}, Iters{"iters", 8},
      VoxelDepth{"voxelDepth", -1}, FullDepth{"fullDepth", DEFAULT_FULL_DEPTH},
#ifndef FAST_COMPILE
      BType{"bType", BOUNDARY_NEUMANN + 1},
#endif  // !FAST_COMPILE
      MaxSolveDepth{"maxSolveDepth"}, Threads{"threads", omp_get_num_procs()};

  cmdLineFloat Color{"color", 16.f}, SamplesPerNode{"samplesPerNode", 1.5f},
      Scale{"scale", 1.1f}, CGSolverAccuracy{"cgAccuracy", float(1e-3)},
      LowResIterMultiplier{"iterMultiplier", 1.f},
      PointWeight{"pointWeight", 4.f};

  cmdLineReadable* params[34] = {
#ifndef FAST_COMPILE
      &Degree,
      &Double,
      &BType,
#endif  // !FAST_COMPILE
      &In,
      &Depth,
      &Out,
      &XForm,
      &Scale,
      &Verbose,
      &CGSolverAccuracy,
      &NoComments,
      &LowResIterMultiplier,
      &KernelDepth,
      &SamplesPerNode,
      &Confidence,
      &NormalWeights,
      &NonManifold,
      &PolygonMesh,
      &ASCII,
      &ShowResidual,
      &VoxelDepth,
      &PointWeight,
      &VoxelGrid,
      &Threads,
      &MaxSolveDepth,
      &AdaptiveExponent,
      &Density,
      &FullDepth,
      &CGDepth,
      &Iters,
      &Color,
      &LinearFit,
      &PrimalVoxel,
      &TempDir,
#if defined(_WIN32) || defined(_WIN64)
      &Performance,
#endif  // _WIN32 || _WIN64
  };
  //  cmdLineReadable* params[34];
  void ShowUsage(char* ex);

  template <class Real>
  struct ColorInfo
  {
    static Point3D<Real> ReadASCII(FILE* fp)
    {
      Point3D<unsigned char> c;
      if (fscanf(fp, " %c %c %c ", &c[0], &c[1], &c[2]) != 3)
        fprintf(stderr, "[ERROR] Failed to read color\n"), exit(0);
      return Point3D<Real>((Real)c[0], (Real)c[1], (Real)c[2]);
    };
    static bool ValidPlyProperties(const bool* props)
    {
      return (props[0] || props[3]) && (props[1] || props[4]) &&
             (props[2] || props[5]);
    }
    const static PlyProperty PlyProperties[];
  };

  double Weight(double v, double start, double end);

#if defined(_WIN32) || defined(_WIN64)
  double PeakMemoryUsageMB(void);
#endif  // _WIN32 || _WIN64

  template <class Real>
  struct OctreeProfiler
  {
    Octree<Real>& tree;
    double t;
    PoissonRecon* Poisson = new PoissonRecon();

    OctreeProfiler(Octree<Real>& t) : tree(t)
    {
      ;
    }
    void start(void)
    {
      t = Time(), tree.resetLocalMemoryUsage();
    }
    void print(const char* header) const
    {
      tree.memoryUsage();
#if defined(_WIN32) || defined(_WIN64)
      if (header)
        printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
               header,
               Time() - t,
               tree.localMemoryUsage(),
               tree.maxMemoryUsage(),
               PeakMemoryUsageMB());
      else
        printf("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
               Time() - t,
               tree.localMemoryUsage(),
               tree.maxMemoryUsage(),
               PeakMemoryUsageMB());
#else   // !_WIN32 && !_WIN64
      if (header)
        printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n",
               header,
               Time() - t,
               tree.localMemoryUsage(),
               tree.maxMemoryUsage());
      else
        printf("%9.1f (s), %9.1f (MB) / %9.1f (MB)\n",
               Time() - t,
               tree.localMemoryUsage(),
               tree.maxMemoryUsage());
#endif  // _WIN32 || _WIN64
    }
    void dumpOutput(const char* header) const
    {
      tree.memoryUsage();
#if defined(_WIN32) || defined(_WIN64)
      if (header)
        Poisson->DumpOutput(
            "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            header,
            Time() - t,
            tree.localMemoryUsage(),
            tree.maxMemoryUsage(),
            PeakMemoryUsageMB());
      else
        Poisson->DumpOutput("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
                            Time() - t,
                            tree.localMemoryUsage(),
                            tree.maxMemoryUsage(),
                            PeakMemoryUsageMB());
#else   // !_WIN32 && !_WIN64
      if (header)
        Poisson->DumpOutput(
            "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            header,
            Time() - t,
            tree.localMemoryUsage(),
            tree.maxMemoryUsage());
      else
        Poisson->DumpOutput("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
                            Time() - t,
                            tree.localMemoryUsage(),
                            tree.maxMemoryUsage());
#endif  // _WIN32 || _WIN64
    }
    void dumpOutput2(std::vector<char*>& comments, const char* header) const
    {
      tree.memoryUsage();
#if defined(_WIN32) || defined(_WIN64)
      if (header)
        Poisson->DumpOutput2(
            comments,
            "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            header,
            Time() - t,
            tree.localMemoryUsage(),
            tree.maxMemoryUsage(),
            PeakMemoryUsageMB());
      else
        Poisson->DumpOutput2(
            comments,
            "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            Time() - t,
            tree.localMemoryUsage(),
            tree.maxMemoryUsage(),
            PeakMemoryUsageMB());
#else   // !_WIN32 && !_WIN64
      if (header)
        Poisson->DumpOutput2(comments,
                             "%s %9.1f (s), %9.1f (MB) / %9.1f (MB)\n",
                             header,
                             Time() - t,
                             tree.localMemoryUsage(),
                             tree.maxMemoryUsage());
      else
        Poisson->DumpOutput2(comments,
                             "%9.1f (s), %9.1f (MB) / %9.1f (MB)\n",
                             Time() - t,
                             tree.localMemoryUsage(),
                             tree.maxMemoryUsage());
#endif  // _WIN32 || _WIN64
    }
  };

  template <class Real>
  XForm4x4<Real> GetPointXForm(OrientedPointStream<Real>& stream,
                               Real scaleFactor)
  {
    Point3D<Real> min, max;
    stream.boundingBox(min, max);
    Point3D<Real> center = (max + min) / 2;
    Real scale           = std::max<Real>(
        max[0] - min[0], std::max<Real>(max[1] - min[1], max[2] - min[2]));
    scale *= scaleFactor;
    for (int i            = 0; i < 3; i++) center[i] -= scale / 2;
    XForm4x4<Real> tXForm = XForm4x4<Real>::Identity(),
                   sXForm = XForm4x4<Real>::Identity();
    for (int i = 0; i < 3; i++)
      sXForm(i, i) = (Real)(1. / scale), tXForm(3, i) = -center[i];
    return sXForm * tXForm;
  }

  template <class Real, int Degree, BoundaryType BType, class Vertex>
  int _Execute(int num_options, char* options[])
  {
    typedef typename Octree<Real>::template DensityEstimator<WEIGHT_DEGREE>
        DensityEstimator;
    typedef typename Octree<Real>::template InterpolationInfo<false>
        InterpolationInfo;
    typedef OrientedPointStream<Real> PointStream;
    typedef OrientedPointStreamWithData<Real, Point3D<Real> >
        PointStreamWithData;
    typedef TransformedOrientedPointStream<Real> XPointStream;
    typedef TransformedOrientedPointStreamWithData<Real, Point3D<Real> >
        XPointStreamWithData;
    Reset<Real>();
    int paramNum = sizeof(params) / sizeof(cmdLineReadable*);
    std::vector<char*> comments;

    if (Verbose.set) echoStdout = 1;

    XForm4x4<Real> xForm, iXForm;
    if (XForm.set)
    {
      FILE* fp = fopen(XForm.value, "r");
      if (!fp)
      {
        fprintf(
            stderr, "[WARNING] Could not read x-form from: %s\n", XForm.value);
        xForm = XForm4x4<Real>::Identity();
      }
      else
      {
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
          {
            float f;
            if (fscanf(fp, " %f ", &f) != 1)
              fprintf(stderr, "[ERROR] Execute: Failed to read xform\n"),
                  exit(0);
            xForm(i, j) = (Real)f;
          }
        fclose(fp);
      }
    }
    else
      xForm = XForm4x4<Real>::Identity();

    DumpOutput2(comments,
                "Running Screened Poisson Reconstruction (Version 9.011)\n");
    char str[1024];
    for (int i = 0; i < paramNum; i++)
    {
      // std::cout << params[i] << std::endl;
      if (params[i]->set)
      {
        params[i]->writeValue(str);
        if (strlen(str))
          DumpOutput2(comments, "\t--%s %s\n", params[i]->name, str);
        else
          DumpOutput2(comments, "\t--%s\n", params[i]->name);
      }
    }
    double startTime = Time();
    Real isoValue    = 0;
    OctNode<TreeNodeData>::SetAllocator(MEMORY_ALLOCATOR_BLOCK_SIZE);
    Octree<Real> tree;
    OctreeProfiler<Real> profiler(tree);
    tree.threads = Threads.value;
    if (!In.set)
    {
      ShowUsage(options[0]);
      return 0;
    }
    if (!MaxSolveDepth.set) MaxSolveDepth.value = Depth.value;

    int kernelDepth = KernelDepth.set ? KernelDepth.value : Depth.value - 2;
    if (kernelDepth > Depth.value)
    {
      fprintf(stderr,
              "[WARNING] %s can't be greater than %s: %d <= %d\n",
              KernelDepth.name,
              Depth.name,
              KernelDepth.value,
              Depth.value);
      kernelDepth = Depth.value;
    }

    int pointCount;

    Real pointWeightSum;
    std::vector<typename Octree<Real>::PointSample>* samples =
        new std::vector<typename Octree<Real>::PointSample>();
    std::vector<ProjectiveData<Point3D<Real>, Real> >* sampleData = NULL;
    DensityEstimator* density = NULL;
    SparseNodeData<Point3D<Real>, NORMAL_DEGREE>* normalInfo = NULL;
    Real targetValue = (Real)0.5;
    // Read in the samples (and color data)
    {
      profiler.start();
      PointStream* pointStream;
      char* ext = GetFileExtension(In.value);
      if (Color.set && Color.value > 0)
      {
        sampleData = new std::vector<ProjectiveData<Point3D<Real>, Real> >();
        if (!strcasecmp(ext, "bnpts"))
          pointStream =
              new BinaryOrientedPointStreamWithData<Real,
                                                    Point3D<Real>,
                                                    float,
                                                    Point3D<unsigned char> >(
                  In.value);
        else if (!strcasecmp(ext, "ply"))
          pointStream =
              new PLYOrientedPointStreamWithData<Real, Point3D<Real> >(
                  In.value,
                  ColorInfo<Real>::PlyProperties,
                  6,
                  ColorInfo<Real>::ValidPlyProperties);
        else
          pointStream =
              new ASCIIOrientedPointStreamWithData<Real, Point3D<Real> >(
                  In.value, ColorInfo<Real>::ReadASCII);
      }
      else
      {
        if (!strcasecmp(ext, "bnpts"))
          pointStream = new BinaryOrientedPointStream<Real, float>(In.value);
        else if (!strcasecmp(ext, "ply"))
          pointStream = new PLYOrientedPointStream<Real>(In.value);
        else
          pointStream = new ASCIIOrientedPointStream<Real>(In.value);
      }
      delete[] ext;
      XPointStream _pointStream(xForm, *pointStream);
      xForm = GetPointXForm(_pointStream, (Real)Scale.value) * xForm;
      if (sampleData)
      {
        XPointStreamWithData _pointStream(xForm,
                                          (PointStreamWithData&)*pointStream);
        pointCount = tree.template init<Point3D<Real> >(
            _pointStream, Depth.value, Confidence.set, *samples, sampleData);
      }
      else
      {
        XPointStream _pointStream(xForm, *pointStream);
        pointCount = tree.template init<Point3D<Real> >(
            _pointStream, Depth.value, Confidence.set, *samples, sampleData);
      }
      iXForm = xForm.inverse();
      delete pointStream;
#pragma omp parallel for num_threads(Threads.value)
      for (int i = 0; i < (int)samples->size(); i++)
        (*samples)[i].sample.data.n *= (Real)-1;

      DumpOutput(
          "Input Points / Samples: %d / %d\n", pointCount, samples->size());
      profiler.dumpOutput2(comments, "# Read input into tree:");
    }
    DenseNodeData<Real, Degree> solution;

    {
      DenseNodeData<Real, Degree> constraints;
      InterpolationInfo* iInfo = NULL;
      int solveDepth           = MaxSolveDepth.value;

      tree.resetNodeIndices();

      // Get the kernel density estimator [If discarding, compute anew.
      // Otherwise,
      // compute once.]
      {
        profiler.start();
        density = tree.template setDensityEstimator<WEIGHT_DEGREE>(
            *samples, kernelDepth, SamplesPerNode.value);
        profiler.dumpOutput2(comments, "#   Got kernel density:");
      }

      // Transform the Hermite samples into a vector field [If discarding,
      // compute
      // anew. Otherwise, compute once.]
      {
        profiler.start();
        normalInfo  = new SparseNodeData<Point3D<Real>, NORMAL_DEGREE>();
        *normalInfo = tree.template setNormalField<NORMAL_DEGREE>(
            *samples, *density, pointWeightSum, BType == BOUNDARY_NEUMANN);
        profiler.dumpOutput2(comments, "#     Got normal field:");
      }

      if (!Density.set) delete density, density = NULL;

      // Trim the tree and prepare for multigrid
      {
        profiler.start();
        std::vector<int> indexMap;

        constexpr int MAX_DEGREE =
            NORMAL_DEGREE > Degree ? NORMAL_DEGREE : Degree;
        tree.template inalizeForBroodedMultigrid<MAX_DEGREE, Degree, BType>(
            FullDepth.value,
            typename Octree<Real>::template HasNormalDataFunctor<NORMAL_DEGREE>(
                *normalInfo),
            &indexMap);

        if (normalInfo) normalInfo->remapIndices(indexMap);
        if (density) density->remapIndices(indexMap);
        profiler.dumpOutput2(comments, "#       Finalized tree:");
      }

      // Add the FEM constraints
      {
        profiler.start();
        constraints = tree.template initDenseNodeData<Degree>();
        tree.template addFEMConstraints<Degree, BType, NORMAL_DEGREE, BType>(
            FEMVFConstraintFunctor<NORMAL_DEGREE, BType, Degree, BType>(1., 0.),
            *normalInfo,
            constraints,
            solveDepth);
        profiler.dumpOutput2(comments, "#  Set FEM constraints:");
      }

      // Free up the normal info [If we don't need it for subseequent
      // iterations.]
      delete normalInfo, normalInfo = NULL;

      // Add the interpolation constraints
      if (PointWeight.value > 0)
      {
        profiler.start();
        iInfo = new InterpolationInfo(tree,
                                      *samples,
                                      targetValue,
                                      AdaptiveExponent.value,
                                      (Real)PointWeight.value * pointWeightSum,
                                      (Real)0);
        tree.template addInterpolationConstraints<Degree, BType>(
            *iInfo, constraints, solveDepth);
        profiler.dumpOutput2(comments, "#Set point constraints:");
      }

      DumpOutput("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n",
                 (int)tree.leaves(),
                 (int)tree.nodes(),
                 (int)tree.ghostNodes());
      DumpOutput("Memory Usage: %.3f MB\n",
                 float(MemoryInfo::Usage()) / (1 << 20));

      // Solve the linear system
      {
        profiler.start();
        typename Octree<Real>::SolverInfo solverInfo;
        solverInfo.cgDepth = CGDepth.value, solverInfo.iters = Iters.value,
        solverInfo.cgAccuracy   = CGSolverAccuracy.value,
        solverInfo.verbose      = Verbose.set,
        solverInfo.showResidual = ShowResidual.set,
        solverInfo.lowResIterMultiplier =
            std::max<double>(1., LowResIterMultiplier.value);
        solution = tree.template solveSystem<Degree, BType>(
            FEMSystemFunctor<Degree, BType>(0, 1., 0),
            iInfo,
            constraints,
            solveDepth,
            solverInfo);
        profiler.dumpOutput2(comments, "# Linear system solved:");
        if (iInfo) delete iInfo, iInfo = NULL;
      }
    }

    char tempHeader[1024];
    {
#if defined(_WIN32) || defined(_WIN64)
      const char FileSeparator = '\\';
#else   // !_WIN
      const char FileSeparator = '/';
#endif  // _WIN
      char tempPath[1024];
      tempPath[0] = 0;
      if (TempDir.set)
        strcpy(tempPath, TempDir.value);
      else
      {
#if defined(_WIN32) || defined(_WIN64)
        GetTempPath(sizeof(tempPath), tempPath);
#else   // !_WIN
#endif  // _WIN
      }
      if (strlen(tempPath) == 0) sprintf(tempPath, ".%c", FileSeparator);
      if (tempPath[strlen(tempPath) - 1] == FileSeparator)
        sprintf(tempHeader, "%sPR_", tempPath);
      else
        sprintf(tempHeader, "%s%cPR_", tempPath, FileSeparator);
    }
    CoredFileMeshData<Vertex> mesh(tempHeader);

    {
      profiler.start();
      double valueSum = 0, weightSum = 0;
      typename Octree<Real>::template MultiThreadedEvaluator<Degree, BType>
          evaluator(&tree, solution, Threads.value);
#pragma omp parallel for num_threads(Threads.value) reduction(+ : valueSum, \
                                                              weightSum)
      for (int j = 0; j < samples->size(); j++)
      {
        ProjectiveData<OrientedPoint3D<Real>, Real>& sample =
            (*samples)[j].sample;
        Real w = sample.weight;
        if (w > 0)
          weightSum += w,
              valueSum += evaluator.value(sample.data.p / sample.weight,
                                          omp_get_thread_num(),
                                          (*samples)[j].node) *
                          w;
      }
      isoValue = (Real)(valueSum / weightSum);
      if (!(Color.set && Color.value > 0) && samples)
        delete samples, samples = NULL;
      profiler.dumpOutput("Got average:");
      DumpOutput("Iso-Value: %e\n", isoValue);
    }

    if (VoxelGrid.set)
    {
      profiler.start();
      FILE* fp = fopen(VoxelGrid.value, "wb");
      if (!fp)
        fprintf(stderr,
                "Failed to open voxel file for writing: %s\n",
                VoxelGrid.value);
      else
      {
        int res              = 0;
        Pointer(Real) values = tree.template voxelEvaluate<Real, Degree, BType>(
            solution, res, isoValue, VoxelDepth.value, PrimalVoxel.set);
        fwrite(&res, sizeof(int), 1, fp);
        if (sizeof(Real) == sizeof(float))
          fwrite(values, sizeof(float), res * res * res, fp);
        else
        {
          float* fValues = new float[res * res * res];
          for (int i   = 0; i < res * res * res; i++)
            fValues[i] = float(values[i]);
          fwrite(fValues, sizeof(float), res * res * res, fp);
          delete[] fValues;
        }
        fclose(fp);
        DeletePointer(values);
      }
      profiler.dumpOutput("Got voxel grid:");
    }

    if (Out.set)
    {
      profiler.start();
      SparseNodeData<ProjectiveData<Point3D<Real>, Real>, DATA_DEGREE>*
          colorData = NULL;
      if (sampleData)
      {
        colorData = new SparseNodeData<ProjectiveData<Point3D<Real>, Real>,
                                       DATA_DEGREE>();
        *colorData = tree.template setDataField<DATA_DEGREE, false>(
            *samples, *sampleData, (DensityEstimator*)NULL);
        delete sampleData, sampleData = NULL;
        for (const OctNode<TreeNodeData>* n = tree.tree().nextNode(); n;
             n                              = tree.tree().nextNode(n))
        {
          ProjectiveData<Point3D<Real>, Real>* clr = (*colorData)(n);
          if (clr) (*clr) *= (Real)pow(Color.value, tree.depth(n));
        }
      }
      tree.template getMCIsoSurface<Degree, BType, WEIGHT_DEGREE, DATA_DEGREE>(
          density,
          colorData,
          solution,
          isoValue,
          mesh,
          !LinearFit.set,
          !NonManifold.set,
          PolygonMesh.set);
      DumpOutput("Vertices / Polygons: %d / %d\n",
                 mesh.outOfCorePointCount() + mesh.inCorePoints.size(),
                 mesh.polygonCount());
      if (PolygonMesh.set)
        profiler.dumpOutput2(comments, "#         Got polygons:");
      else
        profiler.dumpOutput2(comments, "#        Got triangles:");

      if (colorData) delete colorData, colorData = NULL;

      if (NoComments.set)
      {
        if (ASCII.set)
          PlyWritePolygons(Out.value, &mesh, PLY_ASCII, NULL, 0, iXForm);
        else
          PlyWritePolygons(
              Out.value, &mesh, PLY_BINARY_NATIVE, NULL, 0, iXForm);
      }
      else
      {
        if (ASCII.set)
          PlyWritePolygons(Out.value,
                           &mesh,
                           PLY_ASCII,
                           &comments[0],
                           (int)comments.size(),
                           iXForm);
        else
          PlyWritePolygons(Out.value,
                           &mesh,
                           PLY_BINARY_NATIVE,
                           &comments[0],
                           (int)comments.size(),
                           iXForm);
      }
    }
    if (density) delete density, density = NULL;
    DumpOutput2(comments,
                "#          Total Solve: %9.1f (s), %9.1f (MB)\n",
                Time() - startTime,
                tree.maxMemoryUsage());

    return 1;
  }

#if defined(_WIN32) || defined(_WIN64)
  inline double to_seconds(const FILETIME& ft)
  {
    const double low_to_sec  = 100e-9;  // 100 nanoseconds
    const double high_to_sec = low_to_sec * 4294967296.0;
    return ft.dwLowDateTime * low_to_sec + ft.dwHighDateTime * high_to_sec;
  }
#endif  // _WIN32 || _WIN64

#ifndef FAST_COMPILE
  template <class Real, class Vertex>
  int Execute(int num_options, char* options[])
  {
    switch (BType.value)
    {
      case BOUNDARY_FREE + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return _Execute<Real, 1, BOUNDARY_FREE, Vertex>(num_options,
                                                            options);
          case 2:
            return _Execute<Real, 2, BOUNDARY_FREE, Vertex>(num_options,
                                                            options);
          case 3:
            return _Execute<Real, 3, BOUNDARY_FREE, Vertex>(num_options,
                                                            options);
          case 4:
            return _Execute<Real, 4, BOUNDARY_FREE, Vertex>(num_options,
                                                            options);
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 4 are supported");
            return EXIT_FAILURE;
        }
      }
      case BOUNDARY_NEUMANN + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return _Execute<Real, 1, BOUNDARY_NEUMANN, Vertex>(num_options,
                                                               options);
          case 2:
            return _Execute<Real, 2, BOUNDARY_NEUMANN, Vertex>(num_options,
                                                               options);
          case 3:
            return _Execute<Real, 3, BOUNDARY_NEUMANN, Vertex>(num_options,
                                                               options);
          case 4:
            return _Execute<Real, 4, BOUNDARY_NEUMANN, Vertex>(num_options,
                                                               options);
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 4 are supported");
            return EXIT_FAILURE;
        }
      }
      case BOUNDARY_DIRICHLET + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return _Execute<Real, 1, BOUNDARY_DIRICHLET, Vertex>(num_options,
                                                                 options);
          case 2:
            return _Execute<Real, 2, BOUNDARY_DIRICHLET, Vertex>(num_options,
                                                                 options);
          case 3:
            return _Execute<Real, 3, BOUNDARY_DIRICHLET, Vertex>(num_options,
                                                                 options);
          case 4:
            return _Execute<Real, 4, BOUNDARY_DIRICHLET, Vertex>(num_options,
                                                                 options);
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 4 are supported");
            return EXIT_FAILURE;
        }
      }
      default:
        fprintf(stderr, "[ERROR] Not a valid boundary type: %d\n", BType.value);
        return EXIT_FAILURE;
    }
  }
#endif  // !FAST_COMPILE

  int reconstruct(int num_options, char* options[]);
};

template <>
const PlyProperty PoissonRecon::ColorInfo<float>::PlyProperties[] = {
    {const_cast<char*>("r"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[0])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("g"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[1])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("b"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[2])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("red"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[0])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("green"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[1])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("blue"),
     PLY_UCHAR,
     PLY_FLOAT,
     int(offsetof(Point3D<float>, coords[2])),
     0,
     0,
     0,
     0}};
template <>
const PlyProperty PoissonRecon::ColorInfo<double>::PlyProperties[] = {
    {const_cast<char*>("r"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[0])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("g"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[1])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("b"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[2])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("red"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[0])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("green"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[1])),
     0,
     0,
     0,
     0},
    {const_cast<char*>("blue"),
     PLY_UCHAR,
     PLY_DOUBLE,
     int(offsetof(Point3D<double>, coords[2])),
     0,
     0,
     0,
     0}};

#endif  // POISSONRECON_INCLUDED
