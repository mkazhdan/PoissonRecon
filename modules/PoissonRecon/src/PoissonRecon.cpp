#include "Mesh/PoissonRecon/PoissonRecon.h"

void PoissonRecon::DumpOutput(const char* format, ...)
{
  if (outputFile)
  {
    FILE* fp = fopen(outputFile, "a");
    va_list args;
    va_start(args, format);
    vfprintf(fp, format, args);
    fclose(fp);
    va_end(args);
  }
  if (echoStdout)
  {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
  }
}

void PoissonRecon::DumpOutput2(std::vector<char*>& comments,
                               const char* format,
                               ...)
{
  if (outputFile)
  {
    FILE* fp = fopen(outputFile, "a");
    va_list args;
    va_start(args, format);
    vfprintf(fp, format, args);
    fclose(fp);
    va_end(args);
  }

  if (echoStdout)
  {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
  }
  comments.push_back(new char[1024]);
  char* str = comments.back();
  va_list args;
  va_start(args, format);
  vsprintf(str, format, args);
  va_end(args);
  if (str[strlen(str) - 1] == '\n') str[strlen(str) - 1] = 0;
}

void PoissonRecon::ShowUsage(char* ex)
{
  printf("Usage: %s\n", ex);
  printf("\t --%s <input points>\n", In.name);

  printf("\t[--%s <ouput triangle mesh>]\n", Out.name);

  printf("\t[--%s <ouput voxel grid>]\n", VoxelGrid.name);

#ifndef FAST_COMPILE
  printf("\t[--%s <b-spline degree>=%d]\n", Degree.name, Degree.value);

  printf("\t[--%s <boundary type>=%d]\n", BType.name, BType.value);
  for (int i = 0; i < BOUNDARY_COUNT; i++)
    printf("\t\t%d] %s\n", i + 1, BoundaryNames[i]);
#endif  // FAST_COMPILE

  printf(
      "\t[--%s <maximum reconstruction depth>=%d]\n", Depth.name, Depth.value);

  printf("\t[--%s <scale factor>=%f]\n", Scale.name, Scale.value);

  printf("\t[--%s <minimum number of samples per node>=%f]\n",
         SamplesPerNode.name,
         SamplesPerNode.value);

  printf("\t[--%s <interpolation weight>=%.3e]\n",
         PointWeight.name,
         PointWeight.value);

  printf("\t[--%s]\n", Confidence.name);

  printf("\t[--%s]\n", NormalWeights.name);

#ifndef FOR_RELEASE
  printf("\t[--%s <adaptive weighting exponent>=%d]\n",
         AdaptiveExponent.name,
         AdaptiveExponent.value);
#endif  // !FOR_RELEASE

  printf("\t[--%s <iterations>=%d]\n", Iters.name, Iters.value);

#ifndef FOR_RELEASE
  printf("\t[--%s <low-resolution iteration multiplier>=%f]\n",
         LowResIterMultiplier.name,
         LowResIterMultiplier.value);
#endif  // FOR_RELEASE

  printf(
      "\t[--%s <conjugate-gradients depth>=%d]\n", CGDepth.name, CGDepth.value);

#ifndef FOR_RELEASE
  printf("\t[--%s <conjugate-gradients solver accuracy>=%g]\n",
         CGSolverAccuracy.name,
         CGSolverAccuracy.value);
#endif  // !FOR_RELEASE

  printf("\t[--%s <full depth>=%d]\n", FullDepth.name, FullDepth.value);

  printf("\t[--%s <depth at which to extract the voxel grid>=<%s>]\n",
         VoxelDepth.name,
         Depth.name);

  printf("\t[--%s]\n", PrimalVoxel.name);

  printf("\t[--%s <pull factor>]\n", Color.name);

  printf("\t[--%s]\n", Density.name);

  printf("\t[--%s]\n", LinearFit.name);

  printf("\t[--%s]\n", PolygonMesh.name);

#ifndef FOR_RELEASE
  printf("\t[--%s]\n", NonManifold.name);
#endif  // !FOR_RELEASE

#ifdef _OPENMP
  printf("\t[--%s <num threads>=%d]\n", Threads.name, Threads.value);
#endif  // _OPENMP

  printf("\t[--%s]\n", TempDir.name);

  printf("\t[--%s]\n", Verbose.name);

#ifndef FOR_RELEASE
#if defined(_WIN32) || defined(_WIN64)
  printf("\t[--%s]\n", Performance.name);
#endif  // _WIN32 || _WIN64
#endif  // !FOR_RELEASE

#ifndef FOR_RELEASE
  printf("\t[--%s]\n", ASCII.name);

  printf("\t[--%s]\n", NoComments.name);

#ifndef FAST_COMPILE
  printf("\t[--%s]\n", Double.name);
#endif  // FAST_COMPILE
#endif  // !FOR_RELEASE
}

double PoissonRecon::Weight(double v, double start, double end)
{
  v = (v - start) / (end - start);
  if (v < 0)
    return 1.;
  else if (v > 1)
    return 0.;
  else
  {
    // P(x) = a x^3 + b x^2 + c x + d
    //    P (0) = 1 , P (1) = 0 , P'(0) = 0 , P'(1) = 0
    // => d = 1 , a + b + c + d = 0 , c = 0 , 3a + 2b + c = 0
    // => c = 0 , d = 1 , a + b = -1 , 3a + 2b = 0
    // => a = 2 , b = -3 , c = 0 , d = 1
    // => P(x) = 2 x^3 - 3 x^2 + 1
    return 2. * v * v * v - 3. * v * v + 1.;
  }
}

#if defined(_WIN32) || defined(_WIN64)
double PoissonRecon::PeakMemoryUsageMB(void)
{
  HANDLE h = GetCurrentProcess();
  PROCESS_MEMORY_COUNTERS pmc;
  return GetProcessMemoryInfo(h, &pmc, sizeof(pmc))
             ? ((double)pmc.PeakWorkingSetSize) / (1 << 20)
             : 0;
}
#endif  // _WIN32 || _WIN64

int PoissonRecon::reconstruct(int num_options, char* options[])
{
  int i = 0;
  while (i < num_options)
  {
    printf("options[%d]=%s\n", i, options[i++]);
  }
  printf("num_options=%d\n", num_options);

#ifdef ARRAY_DEBUG
  fprintf(stderr, "[WARNING] Running in array debugging mode\n");
#endif  // ARRAY_DEBUG
#if defined(WIN32) && defined(MAX_MEMORY_GB)
  if (MAX_MEMORY_GB > 0)
  {
    SIZE_T peakMemory = 1;
    peakMemory <<= 30;
    peakMemory *= MAX_MEMORY_GB;
    printf("Limiting memory usage to %.2f GB\n", float(peakMemory >> 30));
    HANDLE h = CreateJobObject(NULL, NULL);
    AssignProcessToJobObject(h, GetCurrentProcess());

    JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = {0};
    jeli.BasicLimitInformation.LimitFlags     = JOB_OBJECT_LIMIT_JOB_MEMORY;
    jeli.JobMemoryLimit                       = peakMemory;
    if (!SetInformationJobObject(
            h, JobObjectExtendedLimitInformation, &jeli, sizeof(jeli)))
      fprintf(stderr, "Failed to set memory limit\n");
  }
#endif  // defined( WIN32 ) && defined( MAX_MEMORY_GB )
  double t = Time();

  cmdLineParse(num_options - 1,
               &options[0],
               sizeof(params) / sizeof(cmdLineReadable*),
               params,
               1);
#ifdef FAST_COMPILE
  static const int Degree         = 2;
  static const BoundaryType BType = BOUNDARY_NEUMANN;
  fprintf(stderr,
          "[WARNING] Compiling for degree-%d, boundary-%s, single-precision "
          "_only_\n",
          Degree,
          BoundaryNames[BType]);
  if (Density.set)
    if (Color.set && Color.value > 0)
      return _Execute<float, Degree, BType, PlyColorAndValueVertex<float> >(
          num_options, options);
    else
      return _Execute<float, Degree, BType, PlyValueVertex<float> >(num_options,
                                                                    options);
  else if (Color.set && Color.value > 0)
    return _Execute<float, Degree, BType, PlyColorVertex<float> >(num_options,
                                                                  options);
  else
    return _Execute<float, Degree, BType, PlyVertex<float> >(num_options,
                                                             options);
#else   // !FAST_COMPILE
  {
    if (Density.set)
      if (Color.set && Color.value > 0)
        if (Double.set)
          Execute<double, PlyColorAndValueVertex<float> >(num_options, options);
        else
          Execute<float, PlyColorAndValueVertex<float> >(num_options, options);
      else if (Double.set)
        Execute<double, PlyValueVertex<float> >(num_options, options);
      else
        Execute<float, PlyValueVertex<float> >(num_options, options);
    else if (Color.set && Color.value > 0)
      if (Double.set)
        Execute<double, PlyColorVertex<float> >(num_options, options);
      else
        Execute<float, PlyColorVertex<float> >(num_options, options);
    else if (Double.set)
      Execute<double, PlyVertex<float> >(num_options, options);
    else
      Execute<float, PlyVertex<float> >(num_options, options);
  }
#endif  // FAST_COMPILE
#if defined(_WIN32) || defined(_WIN64)
  if (Performance.set)
  {
    HANDLE cur_thread = GetCurrentThread();
    FILETIME tcreat, texit, tkernel, tuser;
    if (GetThreadTimes(cur_thread, &tcreat, &texit, &tkernel, &tuser))
      printf("Time (Wall/User/Kernel): %.2f / %.2f / %.2f\n",
             Time() - t,
             to_seconds(tuser),
             to_seconds(tkernel));
    else
      printf("Time: %.2f\n", Time() - t);
    HANDLE h = GetCurrentProcess();
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(h, &pmc, sizeof(pmc)))
      printf("Peak Memory (MB): %d\n", (int)(pmc.PeakWorkingSetSize >> 20));
  }
#endif  // _WIN32 || _WIN64
  return EXIT_SUCCESS;
}