
/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer. Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.

Neither the name of the Johns Hopkins University nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef POISSONRECON_INCLUDED
#define POISSONRECON_INCLUDED
#undef SHOW_WARNINGS  // Display compilation warnings
#undef USE_DOUBLE     // If enabled, double-precesion is used
#undef FAST_COMPILE   // If enabled, only a single version of the reconstruction
                      // code is compiled
#undef ARRAY_DEBUG    // If enabled, array access is tested for validity
#define DATA_DEGREE \
  0  // The order of the B-Spline used to splat in data for color interpolation
     // This can be changed to zero if more interpolatory performance is
     // desired.
#define WEIGHT_DEGREE \
  2  // The order of the B-Spline used to splat in the weights for density
     // estimation
#define NORMAL_DEGREE \
  2  // The order of the B-Spline used to splat in the normals for constructing
     // the Laplacian constraints
#define DEFAULT_FEM_DEGREE 2  // The default finite-element degree
#define DEFAULT_FEM_BOUNDARY \
  BOUNDARY_NEUMANN   // The default finite-element boundary type
#define DIMENSION 3  // The dimension of the system

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Mesh/PoissonRecon/MyMiscellany.h"

#include "Mesh/PoissonRecon/CmdLineParser.h"

#include "Mesh/PoissonRecon/PPolynomial.h"

#include "Mesh/PoissonRecon/FEMTree.h"

#include "Mesh/PoissonRecon/Ply.h"

#include "Mesh/PoissonRecon/PointStreamData.h"

class PoissonRecon
{
 public:
  MessageWriter messageWriter;

  const float DefaultPointWeightMultiplier = 2.f;

  cmdLineParameter<char*> In{"in"}, Out{"out"}, TempDir{"tempDir"},
      VoxelGrid{"voxel"}, Tree{"tree"}, Transform{"xForm"};

  cmdLineReadable Performance{"performance"}, ShowResidual{"showResidual"},
      NoComments{"noComments"}, PolygonMesh{"polygonMesh"},
      NonManifold{"nonManifold"}, ASCII{"ascii"}, Density{"density"},
      LinearFit{"linearFit"}, PrimalVoxel{"primalVoxel"},
      ExactInterpolation{"exact"}, Normals{"normals"}, Colors{"colors"},
      Roughness{"roughness"}, Metallic{"metallic"}, Verbose{"verbose"};

  cmdLineParameter<int>
#ifndef FAST_COMPILE
      Degree{"degree", DEFAULT_FEM_DEGREE},
#endif  // !FAST_COMPILE
      Depth{"depth", 8}, KernelDepth{"kernelDepth"}, Iters{"iters", 8},
      FullDepth{"fullDepth", 5}, BaseDepth{"baseDepth", 0},
      BaseVCycles{"baseVCycles", 1},
#ifndef FAST_COMPILE
      BType{"bType", DEFAULT_FEM_BOUNDARY + 1},
#endif  // !FAST_COMPILE
      MaxMemoryGB{"maxMemory", 0}, Threads{"threads", omp_get_num_procs()};

  cmdLineParameter<float> DataX{"data", 32.f},
      SamplesPerNode{"samplesPerNode", 1.5f}, Scale{"scale", 1.1f},
      Width{"width", 0.f}, Confidence{"confidence", 0.f},
      ConfidenceBias{"confidenceBias", 0.f},
      CGSolverAccuracy{"cgAccuracy", 1e-3f}, PointWeight{"pointWeight"};

  cmdLineReadable* params[40] = {
#ifndef FAST_COMPILE
      &Degree,
      &BType,
#endif  // !FAST_COMPILE
      &In,
      &Depth,
      &Out,
      &Transform,
      &Width,
      &Scale,
      &Verbose,
      &CGSolverAccuracy,
      &NoComments,
      &KernelDepth,
      &SamplesPerNode,
      &Confidence,
      &NonManifold,
      &PolygonMesh,
      &ASCII,
      &ShowResidual,
      &ConfidenceBias,
      &BaseDepth,
      &BaseVCycles,
      &PointWeight,
      &VoxelGrid,
      &Threads,
      &Tree,
      &Density,
      &FullDepth,
      &Iters,
      &DataX,
      &Colors,
      &Roughness,
      &Metallic,
      &Normals,
      &LinearFit,
      &PrimalVoxel,
      &TempDir,
      &ExactInterpolation,
      &Performance,
      &MaxMemoryGB,
      NULL};

  void ShowUsage(char* ex)
  {
    printf("Usage: %s\n", ex);
    printf("\t --%s <input points>\n", In.name);
    printf("\t[--%s <ouput triangle mesh>]\n", Out.name);
    printf("\t[--%s <ouput voxel grid>]\n", VoxelGrid.name);
    printf("\t[--%s <ouput fem tree>]\n", Tree.name);
#ifndef FAST_COMPILE
    printf("\t[--%s <b-spline degree>=%d]\n", Degree.name, Degree.value);
    printf("\t[--%s <boundary type>=%d]\n", BType.name, BType.value);
    for (int i = 0; i < BOUNDARY_COUNT; i++)
      printf("\t\t%d] %s\n", i + 1, BoundaryNames[i]);
#endif  // !FAST_COMPILE
    printf("\t[--%s <maximum reconstruction depth>=%d]\n",
           Depth.name,
           Depth.value);
    printf("\t[--%s <voxel width>]\n", Width.name);
    printf("\t[--%s <full depth>=%d]\n", FullDepth.name, FullDepth.value);
    printf("\t[--%s <coarse MG solver depth>=%d]\n",
           BaseDepth.name,
           BaseDepth.value);
    printf("\t[--%s <coarse MG solver v-cycles>=%d]\n",
           BaseVCycles.name,
           BaseVCycles.value);
    printf("\t[--%s <scale factor>=%f]\n", Scale.name, Scale.value);
    printf("\t[--%s <minimum number of samples per node>=%f]\n",
           SamplesPerNode.name,
           SamplesPerNode.value);
    printf("\t[--%s <interpolation weight>=%.3e * <b-spline degree>]\n",
           PointWeight.name,
           DefaultPointWeightMultiplier);
    printf("\t[--%s <iterations>=%d]\n", Iters.name, Iters.value);
    printf("\t[--%s]\n", ExactInterpolation.name);
    printf("\t[--%s <pull factor>=%f]\n", DataX.name, DataX.value);
    printf("\t[--%s]\n", Colors.name);
    printf("\t[--%s]\n", Roughness.name);
    printf("\t[--%s]\n", Metallic.name);
    printf("\t[--%s]\n", Normals.name);
#ifdef _OPENMP
    printf("\t[--%s <num threads>=%d]\n", Threads.name, Threads.value);
#endif  // _OPENMP
    printf("\t[--%s <normal confidence exponent>=%f]\n",
           Confidence.name,
           Confidence.value);
    printf("\t[--%s <normal confidence bias exponent>=%f]\n",
           ConfidenceBias.name,
           ConfidenceBias.value);
    printf("\t[--%s]\n", NonManifold.name);
    printf("\t[--%s]\n", PolygonMesh.name);
    printf("\t[--%s <cg solver accuracy>=%g]\n",
           CGSolverAccuracy.name,
           CGSolverAccuracy.value);
    printf("\t[--%s <maximum memory (in GB)>=%d]\n",
           MaxMemoryGB.name,
           MaxMemoryGB.value);
    printf("\t[--%s]\n", Performance.name);
    printf("\t[--%s]\n", Density.name);
    printf("\t[--%s]\n", LinearFit.name);
    printf("\t[--%s]\n", PrimalVoxel.name);
    printf("\t[--%s]\n", ASCII.name);
    printf("\t[--%s]\n", NoComments.name);
    printf("\t[--%s]\n", TempDir.name);
    printf("\t[--%s]\n", Verbose.name);
  }

  double Weight(double v, double start, double end)
  {
    v = (v - start) / (end - start);
    if (v < 0)
      return 1.;
    else if (v > 1)
      return 0.;
    else
    {
      // P(x) = a x^3 + b x^2 + c x + d
      //		P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
      // =>	d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
      // =>	c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
      // =>	a = 2 , b = -3 , c = 0 , d = 1
      // =>	P(x) = 2 x^3 - 3 x^2 + 1
      return 2. * v * v * v - 3. * v * v + 1.;
    }
  }

  template <unsigned int Dim, class Real>
  struct FEMTreeProfiler
  {
    FEMTree<Dim, Real>& tree;
    double t;
    PoissonRecon* poisson = new PoissonRecon();
    FEMTreeProfiler(FEMTree<Dim, Real>& t) : tree(t)
    {
      ;
    }
    void start(void)
    {
      t = Time(), FEMTree<Dim, Real>::ResetLocalMemoryUsage();
    }
    void print(const char* header) const
    {
      FEMTree<Dim, Real>::MemoryUsage();
      if (header)
        printf("%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
               header,
               Time() - t,
               FEMTree<Dim, Real>::LocalMemoryUsage(),
               FEMTree<Dim, Real>::MaxMemoryUsage(),
               MemoryInfo::PeakMemoryUsageMB());
      else
        printf("%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
               Time() - t,
               FEMTree<Dim, Real>::LocalMemoryUsage(),
               FEMTree<Dim, Real>::MaxMemoryUsage(),
               MemoryInfo::PeakMemoryUsageMB());
    }
    void dumpOutput(const char* header) const
    {
      FEMTree<Dim, Real>::MemoryUsage();
      if (header)
        poisson->messageWriter(
            "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            header,
            Time() - t,
            FEMTree<Dim, Real>::LocalMemoryUsage(),
            FEMTree<Dim, Real>::MaxMemoryUsage(),
            MemoryInfo::PeakMemoryUsageMB());
      else
        poisson->messageWriter(
            "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            Time() - t,
            FEMTree<Dim, Real>::LocalMemoryUsage(),
            FEMTree<Dim, Real>::MaxMemoryUsage(),
            MemoryInfo::PeakMemoryUsageMB());
    }
    void dumpOutput2(std::vector<char*>& comments, const char* header) const
    {
      FEMTree<Dim, Real>::MemoryUsage();
      if (header)
        poisson->messageWriter(
            comments,
            "%s %9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            header,
            Time() - t,
            FEMTree<Dim, Real>::LocalMemoryUsage(),
            FEMTree<Dim, Real>::MaxMemoryUsage(),
            MemoryInfo::PeakMemoryUsageMB());
      else
        poisson->messageWriter(
            comments,
            "%9.1f (s), %9.1f (MB) / %9.1f (MB) / %9.1f (MB)\n",
            Time() - t,
            FEMTree<Dim, Real>::LocalMemoryUsage(),
            FEMTree<Dim, Real>::MaxMemoryUsage(),
            MemoryInfo::PeakMemoryUsageMB());
    }
  };

  template <class Real, unsigned int Dim>
  XForm<Real, Dim + 1> GetBoundingBoxXForm(Point<Real, Dim> min,
                                           Point<Real, Dim> max,
                                           Real scaleFactor)
  {
    Point<Real, Dim> center = (max + min) / 2;
    Real scale              = max[0] - min[0];
    for (int d = 1; d < Dim; d++)
      scale = std::max<Real>(scale, max[d] - min[d]);
    scale *= scaleFactor;
    for (int i = 0; i < Dim; i++) center[i] -= scale / 2;
    XForm<Real, Dim + 1> tXForm = XForm<Real, Dim + 1>::Identity(),
                         sXForm = XForm<Real, Dim + 1>::Identity();
    for (int i = 0; i < Dim; i++)
      sXForm(i, i) = (Real)(1. / scale), tXForm(Dim, i) = -center[i];
    return sXForm * tXForm;
  }
  template <class Real, unsigned int Dim>
  XForm<Real, Dim + 1> GetBoundingBoxXForm(Point<Real, Dim> min,
                                           Point<Real, Dim> max,
                                           Real width,
                                           Real scaleFactor,
                                           int& depth)
  {
    // Get the target resolution (along the largest dimension)
    Real resolution = (max[0] - min[0]) / width;
    for (int d = 1; d < Dim; d++)
      resolution = std::max<Real>(resolution, (max[d] - min[d]) / width);
    resolution *= scaleFactor;
    depth = 0;
    while ((1 << depth) < resolution) depth++;

    Point<Real, Dim> center = (max + min) / 2;
    Real scale              = (1 << depth) * width;

    for (int i = 0; i < Dim; i++) center[i] -= scale / 2;
    XForm<Real, Dim + 1> tXForm = XForm<Real, Dim + 1>::Identity(),
                         sXForm = XForm<Real, Dim + 1>::Identity();
    for (int i = 0; i < Dim; i++)
      sXForm(i, i) = (Real)(1. / scale), tXForm(Dim, i) = -center[i];
    return sXForm * tXForm;
  }

  template <class Real, unsigned int Dim>
  XForm<Real, Dim + 1> GetPointXForm(InputPointStream<Real, Dim>& stream,
                                     Real width,
                                     Real scaleFactor,
                                     int& depth)
  {
    Point<Real, Dim> min, max;
    stream.boundingBox(min, max);
    return GetBoundingBoxXForm(min, max, width, scaleFactor, depth);
  }
  template <class Real, unsigned int Dim>
  XForm<Real, Dim + 1> GetPointXForm(InputPointStream<Real, Dim>& stream,
                                     Real scaleFactor)
  {
    Point<Real, Dim> min, max;
    stream.boundingBox(min, max);
    return GetBoundingBoxXForm(min, max, scaleFactor);
  }

  template <unsigned int Dim, typename Real>
  struct ConstraintDual
  {
    Real target, weight;
    ConstraintDual(Real t, Real w) : target(t), weight(w)
    {
    }
    CumulativeDerivativeValues<Real, Dim, 0> operator()(
        const Point<Real, Dim>& p) const
    {
      return CumulativeDerivativeValues<Real, Dim, 0>(target * weight);
    };
  };
  template <unsigned int Dim, typename Real>
  struct SystemDual
  {
    Real weight;
    SystemDual(Real w) : weight(w)
    {
    }
    CumulativeDerivativeValues<Real, Dim, 0> operator()(
        const Point<Real, Dim>& p,
        const CumulativeDerivativeValues<Real, Dim, 0>& dValues) const
    {
      return dValues * weight;
    };
    CumulativeDerivativeValues<double, Dim, 0> operator()(
        const Point<Real, Dim>& p,
        const CumulativeDerivativeValues<double, Dim, 0>& dValues) const
    {
      return dValues * weight;
    };
  };
  template <unsigned int Dim>
  struct SystemDual<Dim, double>
  {
    typedef double Real;
    Real weight;
    SystemDual(Real w) : weight(w)
    {
    }
    CumulativeDerivativeValues<Real, Dim, 0> operator()(
        const Point<Real, Dim>& p,
        const CumulativeDerivativeValues<Real, Dim, 0>& dValues) const
    {
      return dValues * weight;
    };
  };

  template <typename Vertex,
            typename Real,
            unsigned int... FEMSigs,
            typename... SampleData>
  void ExtractMesh(
      UIntPack<FEMSigs...>,
      std::tuple<SampleData...>,
      FEMTree<sizeof...(FEMSigs), Real>& tree,
      const DenseNodeData<Real, UIntPack<FEMSigs...> >& solution,
      Real isoValue,
      const std::vector<
          typename FEMTree<sizeof...(FEMSigs), Real>::PointSample>* samples,
      std::vector<
          MultiPointStreamData<Real,
                               PointStreamNormal<Real, DIMENSION>,
                               MultiPointStreamData<Real, SampleData...> > >*
          sampleData,
      const typename FEMTree<sizeof...(FEMSigs),
                             Real>::template DensityEstimator<WEIGHT_DEGREE>*
          density,
      std::function<void(
          Vertex&,
          Point<Real, DIMENSION>,
          Real,
          MultiPointStreamData<Real,
                               PointStreamNormal<Real, DIMENSION>,
                               MultiPointStreamData<Real, SampleData...> >)>
          SetVertex,
      std::vector<char*> comments,
      XForm<Real, sizeof...(FEMSigs) + 1> iXForm)
  {
    static const int Dim = sizeof...(FEMSigs);
    typedef UIntPack<FEMSigs...> Sigs;
    typedef PointStreamNormal<Real, Dim> NormalPointSampleData;
    typedef MultiPointStreamData<Real, SampleData...> AdditionalPointSampleData;
    typedef MultiPointStreamData<Real,
                                 NormalPointSampleData,
                                 AdditionalPointSampleData>
        TotalPointSampleData;
    static const unsigned int DataSig =
        FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;
    typedef
        typename FEMTree<Dim, Real>::template DensityEstimator<WEIGHT_DEGREE>
            DensityEstimator;

    // FEMTreeProfiler<Dim, Real> profiler(tree);

    char tempHeader[1024];
    {
      char tempPath[1024];
      tempPath[0] = 0;
      if (TempDir.set)
        strcpy(tempPath, TempDir.value);
      else
        SetTempDirectory(tempPath, sizeof(tempPath));
      if (strlen(tempPath) == 0) sprintf(tempPath, ".%c", FileSeparator);
      if (tempPath[strlen(tempPath) - 1] == FileSeparator)
        sprintf(tempHeader, "%sPR_", tempPath);
      else
        sprintf(tempHeader, "%s%cPR_", tempPath, FileSeparator);
    }
    CoredFileMeshData<Vertex> mesh(tempHeader);

    // profiler.start();
    typename IsoSurfaceExtractor<Dim, Real, Vertex>::IsoStats isoStats;
    if (sampleData)
    {
      SparseNodeData<ProjectiveData<TotalPointSampleData, Real>,
                     IsotropicUIntPack<Dim, DataSig> >
          _sampleData = tree.template setDataField<DataSig, false>(
              *samples, *sampleData, (DensityEstimator*)NULL);
      for (const RegularTreeNode<Dim, FEMTreeNodeData>* n =
               tree.tree().nextNode();
           n;
           n = tree.tree().nextNode(n))
      {
        ProjectiveData<TotalPointSampleData, Real>* clr = _sampleData(n);
        if (clr) (*clr) *= (Real)pow(DataX.value, tree.depth(n));
      }
      isoStats = IsoSurfaceExtractor<Dim, Real, Vertex>::template Extract<
          TotalPointSampleData>(Sigs(),
                                UIntPack<WEIGHT_DEGREE>(),
                                UIntPack<DataSig>(),
                                tree,
                                density,
                                &_sampleData,
                                solution,
                                isoValue,
                                mesh,
                                SetVertex,
                                !LinearFit.set,
                                !NonManifold.set,
                                PolygonMesh.set,
                                false);
    }
    else
      isoStats = IsoSurfaceExtractor<Dim, Real, Vertex>::template Extract<
          TotalPointSampleData>(Sigs(),
                                UIntPack<WEIGHT_DEGREE>(),
                                UIntPack<DataSig>(),
                                tree,
                                density,
                                NULL,
                                solution,
                                isoValue,
                                mesh,
                                SetVertex,
                                !LinearFit.set,
                                !NonManifold.set,
                                PolygonMesh.set,
                                false);
    messageWriter("Vertices / Polygons: %d / %d\n",
                  mesh.outOfCorePointCount() + mesh.inCorePoints.size(),
                  mesh.polygonCount());
    messageWriter(
        "Corners / Vertices / Edges / Surface / Set Table / Copy Finer: %.1f / "
        "%.1f / %.1f / %.1f / %.1f / %.1f (s)\n",
        isoStats.cornersTime,
        isoStats.verticesTime,
        isoStats.edgesTime,
        isoStats.surfaceTime,
        isoStats.setTableTime,
        isoStats.copyFinerTime);
    // if (PolygonMesh.set)
    // profiler.dumpOutput2(comments, "#         Got polygons:");
    // else
    // profiler.dumpOutput2(comments, "#        Got triangles:");

    if (NoComments.set)
      PlyWritePolygons<Vertex, Real, Dim>(
          Out.value,
          &mesh,
          ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE,
          NULL,
          0,
          iXForm);
    else
      PlyWritePolygons<Vertex, Real, Dim>(
          Out.value,
          &mesh,
          ASCII.set ? PLY_ASCII : PLY_BINARY_NATIVE,
          &comments[0],
          (int)comments.size(),
          iXForm);
  }

  template <class Real, typename... SampleData, unsigned int... FEMSigs>
  int Execute(int numOptions, char* options[], UIntPack<FEMSigs...>)
  {
    static const int Dim = sizeof...(FEMSigs);
    typedef UIntPack<FEMSigs...> Sigs;
    typedef UIntPack<FEMSignature<FEMSigs>::Degree...> Degrees;
    typedef UIntPack<
        FEMDegreeAndBType<NORMAL_DEGREE,
                          DerivativeBoundary<FEMSignature<FEMSigs>::BType,
                                             1>::BType>::Signature...>
        NormalSigs;
    static const unsigned int DataSig =
        FEMDegreeAndBType<DATA_DEGREE, BOUNDARY_FREE>::Signature;
    typedef
        typename FEMTree<Dim, Real>::template DensityEstimator<WEIGHT_DEGREE>
            DensityEstimator;
    typedef typename FEMTree<Dim, Real>::template InterpolationInfo<Real, 0>
        InterpolationInfo;
    typedef PointStreamNormal<Real, Dim> NormalPointSampleData;
    typedef MultiPointStreamData<Real, SampleData...> AdditionalPointSampleData;
    typedef MultiPointStreamData<Real,
                                 NormalPointSampleData,
                                 AdditionalPointSampleData>
        TotalPointSampleData;
    typedef InputPointStreamWithData<Real, Dim, TotalPointSampleData>
        InputPointStream;
    typedef TransformedInputPointStreamWithData<Real, Dim, TotalPointSampleData>
        XInputPointStream;
    std::vector<char*> comments;
    messageWriter(
        comments,
        "*************************************************************\n");
    messageWriter(
        comments,
        "*************************************************************\n");
    messageWriter(
        comments,
        "** Running Screened Poisson Reconstruction (Version %s) **\n",
        VERSION);
    messageWriter(
        comments,
        "*************************************************************\n");
    messageWriter(
        comments,
        "*************************************************************\n");

    XForm<Real, Dim + 1> xForm, iXForm;
    if (Transform.set)
    {
      FILE* fp = fopen(Transform.value, "r");
      if (!fp)
      {
        fprintf(stderr,
                "[WARNING] Could not read x-form from: %s\n",
                Transform.value);
        xForm = XForm<Real, Dim + 1>::Identity();
      }
      else
      {
        for (int i = 0; i < Dim + 1; i++)
          for (int j = 0; j < Dim + 1; j++)
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
      xForm = XForm<Real, Dim + 1>::Identity();

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

    double startTime = Time();
    Real isoValue    = 0;

    FEMTree<Dim, Real> tree(MEMORY_ALLOCATOR_BLOCK_SIZE);
    // FEMTreeProfiler<Dim, Real> profiler(tree);

    if (Depth.set && Width.value > 0)
    {
      fprintf(stderr,
              "[WARNING] Both --%s and --%s set, ignoring --%s\n",
              Depth.name,
              Width.name,
              Width.name);
      Width.value = 0;
    }

    int pointCount;

    Real pointWeightSum;
    std::vector<typename FEMTree<Dim, Real>::PointSample>* samples =
        new std::vector<typename FEMTree<Dim, Real>::PointSample>();
    std::vector<TotalPointSampleData>* sampleData            = NULL;
    DensityEstimator* density                                = NULL;
    SparseNodeData<Point<Real, Dim>, NormalSigs>* normalInfo = NULL;
    Real targetValue                                         = (Real)0.5;

    // Read in the samples (and color data)
    {
      // profiler.start();
      InputPointStream* pointStream;
      char* ext  = GetFileExtension(In.value);
      sampleData = new std::vector<TotalPointSampleData>();
      if (!strcasecmp(ext, "bnpts"))
        pointStream =
            new BinaryInputPointStreamWithData<Real, Dim, TotalPointSampleData>(
                In.value, TotalPointSampleData::ReadBinary);
      else if (!strcasecmp(ext, "ply"))
        pointStream =
            new PLYInputPointStreamWithData<Real, Dim, TotalPointSampleData>(
                In.value,
                TotalPointSampleData::PlyReadProperties(),
                TotalPointSampleData::PlyReadNum,
                TotalPointSampleData::ValidPlyReadProperties);
      else
        pointStream =
            new ASCIIInputPointStreamWithData<Real, Dim, TotalPointSampleData>(
                In.value, TotalPointSampleData::ReadASCII);
      delete[] ext;
      typename TotalPointSampleData::Transform _xForm(xForm);
      XInputPointStream _pointStream(
          [&](Point<Real, Dim>& p, TotalPointSampleData& d) {
            p = xForm * p, d = _xForm(d);
          },
          *pointStream);
      if (Width.value > 0)
        xForm =
            GetPointXForm<Real, Dim>(_pointStream,
                                     Width.value,
                                     (Real)(Scale.value > 0 ? Scale.value : 1.),
                                     Depth.value) *
            xForm;
      else
        xForm = Scale.value > 0 ? GetPointXForm<Real, Dim>(_pointStream,
                                                           (Real)Scale.value) *
                                      xForm
                                : xForm;
      {
        typename TotalPointSampleData::Transform _xForm(xForm);
        XInputPointStream _pointStream(
            [&](Point<Real, Dim>& p, TotalPointSampleData& d) {
              p = xForm * p, d = _xForm(d);
            },
            *pointStream);
        auto ProcessDataWithConfidence = [&](const Point<Real, Dim>& p,
                                             TotalPointSampleData& d) {
          Real l = (Real)Length(std::get<0>(d.data).data);
          if (!l || l != l) return (Real)-1.;
          return (Real)pow(l, Confidence.value);
        };
        auto ProcessData = [](const Point<Real, Dim>& p,
                              TotalPointSampleData& d) {
          Real l = (Real)Length(std::get<0>(d.data).data);
          if (!l || l != l) return (Real)-1.;
          std::get<0>(d.data).data /= l;
          return (Real)1.;
        };
        if (Confidence.value > 0)
          pointCount = FEMTreeInitializer<Dim, Real>::template Initialize<
              TotalPointSampleData>(tree.spaceRoot(),
                                    _pointStream,
                                    Depth.value,
                                    *samples,
                                    *sampleData,
                                    true,
                                    tree.nodeAllocator,
                                    tree.initializer(),
                                    ProcessDataWithConfidence);
        else
          pointCount = FEMTreeInitializer<Dim, Real>::template Initialize<
              TotalPointSampleData>(tree.spaceRoot(),
                                    _pointStream,
                                    Depth.value,
                                    *samples,
                                    *sampleData,
                                    true,
                                    tree.nodeAllocator,
                                    tree.initializer(),
                                    ProcessData);
      }
      iXForm = xForm.inverse();
      delete pointStream;

      messageWriter(
          "Input Points / Samples: %d / %d\n", pointCount, samples->size());
      // profiler.dumpOutput2(comments, "# Read input into tree:");
    }
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

    DenseNodeData<Real, Sigs> solution;
    {
      DenseNodeData<Real, Sigs> constraints;
      InterpolationInfo* iInfo = NULL;
      int solveDepth           = Depth.value;

      tree.resetNodeIndices();

      // Get the kernel density estimator
      {
        // profiler.start();
        density = tree.template setDensityEstimator<WEIGHT_DEGREE>(
            *samples, kernelDepth, SamplesPerNode.value, 1);
        // profiler.dumpOutput2(comments, "#   Got kernel density:");
      }

      // Transform the Hermite samples into a vector field
      {
        // profiler.start();
        normalInfo = new SparseNodeData<Point<Real, Dim>, NormalSigs>();
        if (ConfidenceBias.value > 0)
          *normalInfo = tree.setNormalField(
              NormalSigs(),
              *samples,
              *sampleData,
              density,
              pointWeightSum,
              [&](Real conf) {
                return (Real)(log(conf) * ConfidenceBias.value /
                              log(1 << (Dim - 1)));
              });
        else
          *normalInfo = tree.setNormalField(
              NormalSigs(), *samples, *sampleData, density, pointWeightSum);
#pragma omp parallel for
        for (int i = 0; i < normalInfo->size(); i++)
          (*normalInfo)[i] *= (Real)-1.;
        // profiler.dumpOutput2(comments, "#     Got normal field:");
        messageWriter("Point weight / Estimated Area: %g / %g\n",
                      pointWeightSum,
                      pointCount * pointWeightSum);
      }

      if (!Density.set) delete density, density = NULL;
      if (DataX.value <= 0 ||
          (!Colors.set && !Normals.set && !Roughness.set && !Metallic.set))
        delete sampleData, sampleData = NULL;

      // Trim the tree and prepare for multigrid
      {
        // profiler.start();
        constexpr int MAX_DEGREE =
            NORMAL_DEGREE > Degrees::Max() ? NORMAL_DEGREE : Degrees::Max();
        tree.template finalizeForMultigrid<MAX_DEGREE>(
            FullDepth.value,
            typename FEMTree<Dim, Real>::template HasNormalDataFunctor<
                NormalSigs>(*normalInfo),
            normalInfo,
            density);
        // profiler.dumpOutput2(comments, "#       Finalized tree:");
      }
      // Add the FEM constraints
      {
        // profiler.start();
        constraints = tree.initDenseNodeData(Sigs());
        typename FEMIntegrator::template Constraint<Sigs,
                                                    IsotropicUIntPack<Dim, 1>,
                                                    NormalSigs,
                                                    IsotropicUIntPack<Dim, 0>,
                                                    Dim>
            F;
        unsigned int derivatives2[Dim];
        for (int d = 0; d < Dim; d++) derivatives2[d] = 0;
        typedef IsotropicUIntPack<Dim, 1> Derivatives1;
        typedef IsotropicUIntPack<Dim, 0> Derivatives2;
        for (int d = 0; d < Dim; d++)
        {
          unsigned int derivatives1[Dim];
          for (int dd = 0; dd < Dim; dd++) derivatives1[dd] = dd == d ? 1 : 0;
          F.weights[d][TensorDerivatives<Derivatives1>::Index(derivatives1)]
                   [TensorDerivatives<Derivatives2>::Index(derivatives2)] = 1;
        }
        tree.addFEMConstraints(F, *normalInfo, constraints, solveDepth);
        // profiler.dumpOutput2(comments, "#  Set FEM constraints:");
      }

      // Free up the normal info
      delete normalInfo, normalInfo = NULL;

      // Add the interpolation constraints
      if (PointWeight.value > 0)
      {
        // profiler.start();
        if (ExactInterpolation.set)
          iInfo = FEMTree<Dim, Real>::
              template InitializeExactPointInterpolationInfo<Real, 0>(
                  tree,
                  *samples,
                  ConstraintDual<Dim, Real>(
                      targetValue, (Real)PointWeight.value * pointWeightSum),
                  SystemDual<Dim, Real>((Real)PointWeight.value *
                                        pointWeightSum),
                  true,
                  false);
        else
          iInfo = FEMTree<Dim, Real>::
              template InitializeApproximatePointInterpolationInfo<Real, 0>(
                  tree,
                  *samples,
                  ConstraintDual<Dim, Real>(
                      targetValue, (Real)PointWeight.value * pointWeightSum),
                  SystemDual<Dim, Real>((Real)PointWeight.value *
                                        pointWeightSum),
                  true,
                  1);
        tree.addInterpolationConstraints(constraints, solveDepth, *iInfo);
        // profiler.dumpOutput2(comments, "#Set point constraints:");
      }

      messageWriter("Leaf Nodes / Active Nodes / Ghost Nodes: %d / %d / %d\n",
                    (int)tree.leaves(),
                    (int)tree.nodes(),
                    (int)tree.ghostNodes());
      messageWriter("Memory Usage: %.3f MB\n",
                    float(MemoryInfo::Usage()) / (1 << 20));

      // Solve the linear system
      {
        // profiler.start();
        typename FEMTree<Dim, Real>::SolverInfo sInfo;
        sInfo.cgDepth = 0, sInfo.cascadic = true, sInfo.vCycles = 1,
        sInfo.iters = Iters.value, sInfo.cgAccuracy = CGSolverAccuracy.value,
        sInfo.verbose = Verbose.set, sInfo.showResidual = ShowResidual.set,
        sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE,
        sInfo.sliceBlockSize     = 1;
        sInfo.baseDepth          = BaseDepth.value,
        sInfo.baseVCycles        = BaseVCycles.value;
        typename FEMIntegrator::template System<Sigs,
                                                IsotropicUIntPack<Dim, 1> >
            F({0., 1.});
        solution =
            tree.solveSystem(Sigs(), F, constraints, solveDepth, sInfo, iInfo);
        // profiler.dumpOutput2(comments, "# Linear system solved:");
        if (iInfo) delete iInfo, iInfo = NULL;
      }
    }

    {
      // profiler.start();
      double valueSum = 0, weightSum = 0;
      typename FEMTree<Dim, Real>::template MultiThreadedEvaluator<Sigs, 0>
          evaluator(&tree, solution);
#pragma omp parallel for reduction(+ : valueSum, weightSum)
      for (int j = 0; j < samples->size(); j++)
      {
        ProjectiveData<Point<Real, Dim>, Real>& sample = (*samples)[j].sample;
        Real w                                         = sample.weight;
        if (w > 0)
          weightSum += w,
              valueSum += evaluator.values(sample.data / sample.weight,
                                           omp_get_thread_num(),
                                           (*samples)[j].node)[0] *
                          w;
      }
      isoValue = (Real)(valueSum / weightSum);
      if (DataX.value <= 0 ||
          (!Colors.set && !Normals.set && !Roughness.set && !Metallic.set))
        delete samples, samples = NULL;
      // profiler.dumpOutput("Got average:");
      messageWriter("Iso-Value: %e = %g / %g\n", isoValue, valueSum, weightSum);
    }
    if (Tree.set)
    {
      FILE* fp = fopen(Tree.value, "wb");
      if (!fp)
        fprintf(stderr,
                "[ERROR] Failed to open file for writing: %s\n",
                Tree.value),
            exit(0);
      FEMTree<Dim, Real>::WriteParameter(fp);
      DenseNodeData<Real, Sigs>::WriteSignatures(fp);
      tree.write(fp);
      solution.write(fp);
      fclose(fp);
    }

    if (VoxelGrid.set)
    {
      FILE* fp = fopen(VoxelGrid.value, "wb");
      if (!fp)
        fprintf(stderr,
                "Failed to open voxel file for writing: %s\n",
                VoxelGrid.value);
      else
      {
        int res = 0;
        // profiler.start();
        Pointer(Real) values = tree.template regularGridEvaluate<true>(
            solution, res, -1, PrimalVoxel.set);
#pragma omp parallel for
        for (int i = 0; i < res * res * res; i++) values[i] -= isoValue;
        // profiler.dumpOutput("Got voxel grid:");
        fwrite(&res, sizeof(int), 1, fp);
        if (typeid(Real) == typeid(float))
          fwrite(values, sizeof(float), res * res * res, fp);
        else
        {
          float* fValues = new float[res * res * res];
          for (int i = 0; i < res * res * res; i++)
            fValues[i] = float(values[i]);
          fwrite(fValues, sizeof(float), res * res * res, fp);
          delete[] fValues;
        }
        fclose(fp);
        DeletePointer(values);
      }
    }

    if (Out.set)
    {
      if (Normals.set)
      {
        if (Density.set)
        {
          typedef PlyVertexWithData<
              Real,
              Dim,
              MultiPointStreamData<Real,
                                   PointStreamNormal<Real, Dim>,
                                   PointStreamValue<Real>,
                                   AdditionalPointSampleData> >
              Vertex;
          std::function<void(
              Vertex&, Point<Real, Dim>, Real, TotalPointSampleData)>
              SetVertex = [](Vertex& v,
                             Point<Real, Dim> p,
                             Real w,
                             TotalPointSampleData d) {
                v.point = p, std::get<0>(v.data.data) = std::get<0>(d.data),
                std::get<1>(v.data.data).data = w,
                std::get<2>(v.data.data)      = std::get<1>(d.data);
              };
          ExtractMesh<Vertex>(UIntPack<FEMSigs...>(),
                              std::tuple<SampleData...>(),
                              tree,
                              solution,
                              isoValue,
                              samples,
                              sampleData,
                              density,
                              SetVertex,
                              comments,
                              iXForm);
        }
        else
        {
          typedef PlyVertexWithData<
              Real,
              Dim,
              MultiPointStreamData<Real,
                                   PointStreamNormal<Real, Dim>,
                                   AdditionalPointSampleData> >
              Vertex;
          std::function<void(
              Vertex&, Point<Real, Dim>, Real, TotalPointSampleData)>
              SetVertex = [](Vertex& v,
                             Point<Real, Dim> p,
                             Real w,
                             TotalPointSampleData d) {
                v.point = p, std::get<0>(v.data.data) = std::get<0>(d.data),
                std::get<1>(v.data.data) = std::get<1>(d.data);
              };
          ExtractMesh<Vertex>(UIntPack<FEMSigs...>(),
                              std::tuple<SampleData...>(),
                              tree,
                              solution,
                              isoValue,
                              samples,
                              sampleData,
                              density,
                              SetVertex,
                              comments,
                              iXForm);
        }
      }
      else
      {
        if (Density.set)
        {
          typedef PlyVertexWithData<
              Real,
              Dim,
              MultiPointStreamData<Real,
                                   PointStreamValue<Real>,
                                   AdditionalPointSampleData> >
              Vertex;
          std::function<void(
              Vertex&, Point<Real, Dim>, Real, TotalPointSampleData)>
              SetVertex = [](Vertex& v,
                             Point<Real, Dim> p,
                             Real w,
                             TotalPointSampleData d) {
                v.point = p, std::get<0>(v.data.data).data = w,
                std::get<1>(v.data.data) = std::get<1>(d.data);
              };
          ExtractMesh<Vertex>(UIntPack<FEMSigs...>(),
                              std::tuple<SampleData...>(),
                              tree,
                              solution,
                              isoValue,
                              samples,
                              sampleData,
                              density,
                              SetVertex,
                              comments,
                              iXForm);
        }
        else
        {
          typedef PlyVertexWithData<
              Real,
              Dim,
              MultiPointStreamData<Real, AdditionalPointSampleData> >
              Vertex;
          std::function<void(
              Vertex&, Point<Real, Dim>, Real, TotalPointSampleData)>
              SetVertex = [](Vertex& v,
                             Point<Real, Dim> p,
                             Real w,
                             TotalPointSampleData d) {
                v.point = p, std::get<0>(v.data.data) = std::get<1>(d.data);
              };
          ExtractMesh<Vertex>(UIntPack<FEMSigs...>(),
                              std::tuple<SampleData...>(),
                              tree,
                              solution,
                              isoValue,
                              samples,
                              sampleData,
                              density,
                              SetVertex,
                              comments,
                              iXForm);
        }
      }
      if (sampleData)
      {
        delete sampleData;
        sampleData = NULL;
      }
    }
    if (density) delete density, density = NULL;
    messageWriter(comments,
                  "#          Total Solve: %9.1f (s), %9.1f (MB)\n",
                  Time() - startTime,
                  FEMTree<Dim, Real>::MaxMemoryUsage());

    return 1;
  }

#ifndef FAST_COMPILE
  template <unsigned int Dim, class Real, typename... SampleData>
  int Execute(int numOptions, char* options[])
  {
    switch (BType.value)
    {
      case BOUNDARY_FREE + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<1, BOUNDARY_FREE>::Signature>());
          case 2:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<2, BOUNDARY_FREE>::Signature>());
            //				case 3: return Execute< Real ,
            // SampleData
            //...
            //>( numOptions , options , IsotropicUIntPack< Dim ,
            // FEMDegreeAndBType< 3 ,
            // BOUNDARY_FREE >::Signature >() ); case 4: return Execute< Real ,
            // SampleData ... >( numOptions , options , IsotropicUIntPack< Dim ,
            // FEMDegreeAndBType< 4 , BOUNDARY_FREE
            // >::Signature >() );
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 2 are supported");
            return EXIT_FAILURE;
        }
      }
      case BOUNDARY_NEUMANN + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature>());
          case 2:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<2, BOUNDARY_NEUMANN>::Signature>());
            //				case 3: return Execute< Real ,
            // SampleData
            //...
            //>( numOptions , options , IsotropicUIntPack< Dim ,
            // FEMDegreeAndBType< 3 ,
            // BOUNDARY_NEUMANN >::Signature >() ); case 4: return Execute< Real
            // , SampleData ... >( numOptions , options , IsotropicUIntPack< Dim
            // , FEMDegreeAndBType< 4 , BOUNDARY_NEUMANN
            // >::Signature >() );
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 2 are supported");
            return EXIT_FAILURE;
        }
      }
      case BOUNDARY_DIRICHLET + 1:
      {
        switch (Degree.value)
        {
          case 1:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<1, BOUNDARY_DIRICHLET>::Signature>());
          case 2:
            return Execute<Real, SampleData...>(
                numOptions,
                options,
                IsotropicUIntPack<
                    Dim,
                    FEMDegreeAndBType<2, BOUNDARY_DIRICHLET>::Signature>());
            //			case 3: return Execute< Real , SampleData ... >(
            // numOptions , options , IsotropicUIntPack< Dim ,
            // FEMDegreeAndBType< 3 , BOUNDARY_DIRICHLET >::Signature >() );
            // case 4: return Execute< Real , SampleData ... >( numOptions ,
            // options , IsotropicUIntPack< Dim , FEMDegreeAndBType< 4 ,
            // BOUNDARY_DIRICHLET >::Signature >() );
          default:
            fprintf(stderr,
                    "[ERROR] Only B-Splines of degree 1 - 2 are supported");
            return EXIT_FAILURE;
        }
      }
      default:
        fprintf(stderr, "[ERROR] Not a valid boundary type: %d\n", BType.value);
        return EXIT_FAILURE;
    }
  }
#endif  // !FAST_COMPILE

  int reconstruct(int numOptions, char* options[])
  {
    Timer timer;
#ifdef ARRAY_DEBUG
    fprintf(stderr, "[WARNING] Array debugging enabled\n");
#endif  // ARRAY_DEBUG

    cmdLineParse(numOptions - 1, &options[1], params);
    if (MaxMemoryGB.value > 0) SetPeakMemoryMB(MaxMemoryGB.value << 10);
    omp_set_num_threads(Threads.value > 1 ? Threads.value : 1);
    messageWriter.echoSTDOUT = Verbose.set;

    if (!In.set)
    {
      ShowUsage(options[0]);
      return 0;
    }
    if (DataX.value <= 0)
      Normals.set = Colors.set = Roughness.set = Metallic.set = false;
    if (BaseDepth.value > FullDepth.value)
    {
      if (BaseDepth.set)
        fprintf(
            stderr,
            "[WARNING] Base depth must be smaller than full depth: %d <= %d\n",
            BaseDepth.value,
            FullDepth.value);
      BaseDepth.value = FullDepth.value;
    }

#ifdef USE_DOUBLE
    typedef double Real;
#else   // !USE_DOUBLE
    typedef float Real;
#endif  // USE_DOUBLE

#ifdef FAST_COMPILE
    static const int Degree         = DEFAULT_FEM_DEGREE;
    static const BoundaryType BType = DEFAULT_FEM_BOUNDARY;
    typedef IsotropicUIntPack<DIMENSION,
                              FEMDegreeAndBType<Degree, BType>::Signature>
        FEMSigs;
    fprintf(
        stderr,
        "[WARNING] Compiled for degree-%d, boundary-%s, %s-precision _only_\n",
        Degree,
        BoundaryNames[BType],
        sizeof(Real) == 4 ? "single" : "double");
    if (!PointWeight.set)
      PointWeight.value = DefaultPointWeightMultiplier * Degree;
    if (Colors.set)
    {
      if (!Roughness.set && !Metallic.set)
        Execute<Real, PointStreamColor<DefaultFloatType> >(
            numOptions, options, FEMSigs());

      if (Roughness.set && !Metallic.set)
        Execute<Real,
                PointStreamColor<DefaultFloatType>,
                PointStreamRoughness<DefaultFloatType> >(
            numOptions, options, FEMSigs());

      if (!Roughness.set && Metallic.set)
        Execute<Real,
                PointStreamColor<DefaultFloatType>,
                PointStreamMetal<DefaultFloatType> >(
            numOptions, options, FEMSigs());

      if (Roughness.set && Metallic.set)
        Execute<Real,
                PointStreamColor<DefaultFloatType>,
                PointStreamRoughness<DefaultFloatType>,
                PointStreamMetal<DefaultFloatType> >(
            numOptions, options, FEMSigs());
    }
    else
      Execute<DefaultFloatType>(numOptions, options, FEMSigs());
#else   // !FAST_COMPILE
    if (!PointWeight.set)
      PointWeight.value = DefaultPointWeightMultiplier * Degree.value;
    if (Colors.set)
    {
      if (!Roughness.set && !Metallic.set)
        Execute<DIMENSION, float, PointStreamColor<float> >(numOptions,
                                                            options);

      if (Roughness.set && !Metallic.set)
        Execute<DIMENSION,
                float,
                PointStreamColor<float>,
                PointStreamRoughness<float> >(numOptions, options);

      if (!Roughness.set && Metallic.set)
        Execute<DIMENSION,
                float,
                PointStreamColor<float>,
                PointStreamMetal<float> >(numOptions, options);

      if (Roughness.set && Metallic.set)
        Execute<DIMENSION,
                float,
                PointStreamColor<float>,
                PointStreamRoughness<float>,
                PointStreamMetal<float> >(numOptions, options);
    }
    else
      Execute<DIMENSION, float>(numOptions, options);
#endif  // FAST_COMPILE
    if (Performance.set)
    {
      printf(
          "Time (Wall/CPU): %.2f / %.2f\n", timer.wallTime(), timer.cpuTime());
      printf("Peak Memory (MB): %d\n", MemoryInfo::PeakMemoryUsageMB());
    }
    return EXIT_SUCCESS;
  }
};

#endif  // POISSONRECON_INCLUDED