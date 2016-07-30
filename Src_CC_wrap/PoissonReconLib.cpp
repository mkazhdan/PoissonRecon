//##########################################################################
//#                                                                        #
//#               CLOUDCOMPARE WRAPPER: PoissonReconLib                    #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 or later of the License.      #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#               COPYRIGHT: Daniel Girardeau-Montaut                      #
//#                                                                        #
//##########################################################################

#include "PoissonReconLib.h"

#ifdef _WIN32
#include <Windows.h>
#include <Psapi.h>
#endif // _WIN32

#ifdef WITH_OPENMP
#include <omp.h>
#endif

//PoissonRecon
#include "../Src/MemoryUsage.h"
#include "../Src/MyTime.h"
#include "../Src/Ply.h"
#include "../Src/Array.h"
#include "../Src/Octree.h"
#include "../Src/SparseMatrix.h"

#define DumpOutput(...) ((void)0)
#include "../Src/MultiGridOctreeData.h" //only after DumpOutput has been defined!

#define BSPLINE_DEGREE 2

PoissonReconLib::Parameters::Parameters()
	: depth(8) //8
	, cgDepth(0) //0
	, kernelDepth(0) //?
	, adaptiveExp(1) //AdaptiveExponent (1)
	, iters(8) //8
	, fullDepth(5) //5
	, minDepth(0) //0
	, maxSolveDepth(0) //?
	, dirichlet(true) //true
	, threads(1) //ideally omp_get_num_procs()
	, samplesPerNode(1.5f) //1.5f
	, scale(1.1f) //1.1f
	, cgAccuracy(1.0e-3f) //1.0e-3f
	, pointWeight(4.0f) //4.0f
	, complete(false)
	, showResidual(false)
	, confidence(false)
	, normalWeights(false)
	, nonManifold(false)
	, density(false)
	, colorInterp(16.0f)
{
#ifdef WITH_OPENMP
	threads = omp_get_num_procs();
#endif
}

template< class PointCoordinateType, class Real, int Degree, class Vertex >
bool Execute(PoissonReconLib::Parameters params, OrientedPointStream< PointCoordinateType >* pointStream, CoredVectorMeshData< Vertex >& mesh)
{
	XForm4x4< Real > xForm = XForm4x4< Real >::Identity();
	XForm4x4< Real > iXForm = xForm.inverse();

	//DGM: reset static parameters!!!
	TreeNodeData::NodeCount = 0;

	Octree< Real > tree;
	tree.threads = params.threads;

	if (params.maxSolveDepth == 0)
		params.maxSolveDepth = params.depth;

	OctNode< TreeNodeData >::SetAllocator( MEMORY_ALLOCATOR_BLOCK_SIZE );

	if (params.maxSolveDepth < 2)
		return false;
	int kernelDepth = params.kernelDepth != 0 ?  params.kernelDepth : params.maxSolveDepth-2;
	if( kernelDepth > params.depth )
		return false;
	params.fullDepth = std::min(params.fullDepth, params.depth);

	tree.maxMemoryUsage = 0;
	SparseNodeData< PointData< Real > , 0 >* pointInfo = new SparseNodeData< PointData < Real > , 0 >();
	SparseNodeData< Point3D< Real > , NORMAL_DEGREE >* normalInfo = new SparseNodeData< Point3D< Real > , NORMAL_DEGREE >();
	SparseNodeData< Real , WEIGHT_DEGREE >* densityWeights = new SparseNodeData< Real , WEIGHT_DEGREE >();
	SparseNodeData< Real , NORMAL_DEGREE >* nodeWeights = new SparseNodeData< Real , NORMAL_DEGREE >();
	typedef typename Octree< Real >::template ProjectiveData< Point3D< Real > > ProjectiveColor;
	SparseNodeData< ProjectiveColor , DATA_DEGREE >* colorData = NULL;

	int pointCount = tree.template SetTree< PointCoordinateType, NORMAL_DEGREE , WEIGHT_DEGREE , DATA_DEGREE , Point3D< unsigned char > >(
									pointStream,
									params.minDepth,
									params.depth,
									params.fullDepth,
									kernelDepth,
									static_cast<Real>(params.samplesPerNode),
									params.scale,
									params.confidence,
									params.normalWeights,
									params.pointWeight,
									params.adaptiveExp,
									*densityWeights,
									*pointInfo,
									*normalInfo,
									*nodeWeights,
									colorData,
									xForm,
									params.dirichlet,
									params.complete );

	if( !params.density )
	{
		delete densityWeights;
		densityWeights = NULL;
	}
	//reamp indexes
	{
		std::vector< int > indexMap;
		if( NORMAL_DEGREE > Degree )
			tree.template EnableMultigrid< NORMAL_DEGREE >( &indexMap );
		else
			tree.template EnableMultigrid<        Degree >( &indexMap );

		if (pointInfo)
			pointInfo->remapIndices( indexMap );
		if (normalInfo)
			normalInfo->remapIndices( indexMap );
		if (densityWeights)
			densityWeights->remapIndices( indexMap );
		if (nodeWeights)
			nodeWeights->remapIndices( indexMap );
	}

	double maxMemoryUsage = tree.maxMemoryUsage;
	tree.maxMemoryUsage = 0;
	DenseNodeData< Real , Degree > constraints = tree.template SetLaplacianConstraints< Degree >( *normalInfo );
	delete normalInfo;
	normalInfo = NULL;
	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );

	tree.maxMemoryUsage = 0;
	DenseNodeData< Real , Degree > solution = tree.SolveSystem( *pointInfo , constraints , params.showResidual , params.iters, params.maxSolveDepth, params.cgDepth, params.cgAccuracy );
	delete pointInfo;
	pointInfo = NULL;
	constraints.resize(0);
	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );

	Real isoValue = tree.GetIsoValue( solution , *nodeWeights );
	delete nodeWeights;
	nodeWeights = NULL;
	//DumpOutput( "Iso-Value: %e\n" , isoValue );

	//output
	tree.maxMemoryUsage = 0;
	tree.template GetMCIsoSurface< Degree , WEIGHT_DEGREE , DATA_DEGREE >(
							densityWeights,
							colorData,
							solution,
							isoValue,
							mesh,
							true,
							!params.nonManifold,
							false );

	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );

	//DumpOutput( "Vertices / Polygons: %d / %d\n" , mesh.outOfCorePointCount()+mesh.inCorePoints.size() , mesh.polygonCount() );

	solution.resize(0);

	return true;
}

template< class PointCoordinateType, class Real, int Degree, class Vertex >
bool Execute(PoissonReconLib::Parameters params, OrientedPointStreamWithData< PointCoordinateType , Point3D< unsigned char > >* pointStream, CoredVectorMeshData< Vertex >& mesh)
{
	XForm4x4< Real > xForm = XForm4x4< Real >::Identity();
	XForm4x4< Real > iXForm = xForm.inverse();

	//DGM: reset static parameters!!!
	TreeNodeData::NodeCount = 0;

	Octree< Real > tree;
	tree.threads = params.threads;

	if (params.maxSolveDepth == 0)
		params.maxSolveDepth = params.depth;

	OctNode< TreeNodeData >::SetAllocator( MEMORY_ALLOCATOR_BLOCK_SIZE );

	if (params.maxSolveDepth < 2)
		return false;
	int kernelDepth = params.kernelDepth != 0 ?  params.kernelDepth : params.maxSolveDepth-2;
	if( kernelDepth > params.depth )
		return false;
	params.fullDepth = std::min(params.fullDepth, params.depth);

	tree.maxMemoryUsage = 0;
	SparseNodeData< PointData< Real > , 0 >* pointInfo = new SparseNodeData< PointData < Real > , 0 >();
	SparseNodeData< Point3D< Real > , NORMAL_DEGREE >* normalInfo = new SparseNodeData< Point3D< Real > , NORMAL_DEGREE >();
	SparseNodeData< Real , WEIGHT_DEGREE >* densityWeights = new SparseNodeData< Real , WEIGHT_DEGREE >();
	SparseNodeData< Real , NORMAL_DEGREE >* nodeWeights = new SparseNodeData< Real , NORMAL_DEGREE >();
	typedef typename Octree< Real >::template ProjectiveData< Point3D< Real > > ProjectiveColor;
	SparseNodeData< ProjectiveColor , DATA_DEGREE > colorData;


	int pointCount = tree.template SetTree< PointCoordinateType, NORMAL_DEGREE , WEIGHT_DEGREE , DATA_DEGREE , Point3D< unsigned char > >
								(
									pointStream,
									params.minDepth,
									params.depth,
									params.fullDepth,
									kernelDepth,
									static_cast<Real>(params.samplesPerNode),
									params.scale,
									params.confidence,
									params.normalWeights,
									params.pointWeight,
									params.adaptiveExp,
									*densityWeights,
									*pointInfo,
									*normalInfo,
									*nodeWeights,
									&colorData,
									xForm,
									params.dirichlet,
									params.complete );

	for (const OctNode< TreeNodeData >* n = tree.tree().nextNode(); n != NULL; n = tree.tree().nextNode( n ) )
	{
		int idx = colorData.index(n);
		if (idx >= 0)
			colorData.data[idx] *= static_cast<Real>(pow(params.colorInterp, n->depth()));
	}

	if( !params.density )
	{
		delete densityWeights;
		densityWeights = NULL;
	}
	//reamp indexes
	{
		std::vector< int > indexMap;
		if( NORMAL_DEGREE > Degree )
			tree.template EnableMultigrid< NORMAL_DEGREE >( &indexMap );
		else
			tree.template EnableMultigrid<        Degree >( &indexMap );

		if (pointInfo)
			pointInfo->remapIndices( indexMap );
		if (normalInfo)
			normalInfo->remapIndices( indexMap );
		if (densityWeights)
			densityWeights->remapIndices( indexMap );
		if (nodeWeights)
			nodeWeights->remapIndices( indexMap );
		colorData.remapIndices( indexMap );
	}

	//DumpOutput( "Input Points: %d\n" , pointCount );
	//DumpOutput( "Leaves/Nodes: %d/%d\n" , tree.tree.leaves() , tree.tree.nodes() );

	double maxMemoryUsage = tree.maxMemoryUsage;
	tree.maxMemoryUsage = 0;

	DenseNodeData< Real , Degree > constraints = tree.template SetLaplacianConstraints< Degree >( *normalInfo );
	delete normalInfo;
	normalInfo = 0;

	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );
	tree.maxMemoryUsage = 0;

	DenseNodeData< Real , Degree > solution = tree.SolveSystem( *pointInfo , constraints , params.showResidual , params.iters, params.maxSolveDepth, params.cgDepth, params.cgAccuracy );

	delete pointInfo;
	pointInfo = 0;
	constraints.resize(0);

	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );

	Real isoValue = tree.GetIsoValue( solution , *nodeWeights );
	delete nodeWeights;
	nodeWeights = 0;

	//DumpOutput( "Iso-Value: %e\n" , isoValue );

	//output
	tree.maxMemoryUsage = 0;

	tree.template GetMCIsoSurface< Degree , WEIGHT_DEGREE , DATA_DEGREE >(
							densityWeights ? GetPointer( *densityWeights ) : NullPointer( Real ),
							&colorData,
							solution,
							isoValue,
							mesh,
							true,
							!params.nonManifold,
							false );

	maxMemoryUsage = std::max< double >( maxMemoryUsage , tree.maxMemoryUsage );

	//DumpOutput( "Vertices / Polygons: %d / %d\n" , mesh.outOfCorePointCount()+mesh.inCorePoints.size() , mesh.polygonCount() );

	solution.resize(0);

	return true;
}

bool PoissonReconLib::Reconstruct(Parameters params, OrientedPointStreamWithData< float , Point3D< unsigned char > >* pointStream, CoredVectorMeshData< PlyColorAndValueVertex< float > >& mesh)
{
	return Execute<	float,
					float,
					BSPLINE_DEGREE,
					PlyColorAndValueVertex< float > > (params, pointStream, mesh);
}

bool PoissonReconLib::Reconstruct(Parameters params, OrientedPointStream< float >* pointStream, CoredVectorMeshData< PlyValueVertex< float > >& mesh)
{
	return Execute<	float,
					float,
					BSPLINE_DEGREE,
					PlyValueVertex< float > > (params, pointStream, mesh);
}

bool PoissonReconLib::Reconstruct(Parameters params, OrientedPointStreamWithData< double , Point3D< unsigned char > >* pointStream, CoredVectorMeshData< PlyColorAndValueVertex< double > >& mesh)
{
	return Execute<	double,
					double,
					BSPLINE_DEGREE,
					PlyColorAndValueVertex< double > > (params, pointStream, mesh);
}

bool PoissonReconLib::Reconstruct(Parameters params, OrientedPointStream< double >* pointStream, CoredVectorMeshData< PlyValueVertex< double > >& mesh)
{
	return Execute<	double,
					double,
					BSPLINE_DEGREE,
					PlyValueVertex< double > > (params, pointStream, mesh);
}
