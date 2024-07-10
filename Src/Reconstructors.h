/*
Copyright (c) 2022, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef RECONSTRUCTORS_INCLUDED
#define RECONSTRUCTORS_INCLUDED

#include "PreProcessor.h"
#include "MyMiscellany.h"
#include "DataStream.imp.h"
#include "FEMTree.h"
#include "PointExtent.h"

namespace Reconstructor
{
	static const unsigned int DataDegree = 0;	// The order of the B-Spline used to splat in data for auxiliary data interpolation
	static const unsigned int WeightDegree = 2;	// The order of the B-Spline used to splat in the weights for density estimation

	// For clarity, to distinguies betwen the case that a Point is referencing a position and a normal
	template< typename Real , unsigned int Dim > using Position = Point< Real , Dim >;
	template< typename Real , unsigned int Dim > using Normal = Point< Real , Dim >;
	template< typename Real , unsigned int Dim > using Gradient = Point< Real , Dim >;

#include "Reconstructors.streams.h"

	// Declare a type for storing the solution information
	template< typename Real , unsigned int Dim , unsigned int FEMSig , typename ... AuxData > struct Implicit;

	// Parameters for mesh extraction
	struct LevelSetExtractionParameters
	{
		bool linearFit;
		bool outputGradients;
		bool forceManifold;
		bool polygonMesh;
		bool verbose;
		LevelSetExtractionParameters( void ) : linearFit(false) , outputGradients(false) , forceManifold(true) , polygonMesh(false) , verbose(false) {}
	};

	// "Private" function for extracting meshes
	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename OutputVertexStream , typename ImplicitType , unsigned int ... FEMSigs >
	void _ExtractLevelSet( UIntPack< FEMSigs ... >  , const ImplicitType &implicit , OutputVertexStream &vertexStream , OutputFaceStream< Dim-1 > &faceStream , typename Implicit< Real , Dim , FEMSig >::LevelSetExtractionParameters params );

	// Specialized solution information without auxiliary data
	template< typename Real , unsigned int Dim , unsigned int FEMSig >
	struct Implicit< Real , Dim , FEMSig >
	{
		// The signature pack
		typedef IsotropicUIntPack< Dim , FEMSig > Sigs;

		// The type representing the point sampling density
		typedef typename FEMTree< Dim , Real >::template DensityEstimator< Reconstructor::WeightDegree > DensityEstimator;

		// The constructor
		Implicit( void ) : density(NULL) , isoValue(0) , tree(MEMORY_ALLOCATOR_BLOCK_SIZE) , unitCubeToModel( XForm< Real , Dim+1 >::Identity() ){}

		// The desctructor
		~Implicit( void ){ delete density ; density = NULL; }

		// The transformation taking points in the unit cube back to world coordinates
		XForm< Real , Dim+1 > unitCubeToModel;

		// The octree adapted to the points
		FEMTree< Dim , Real > tree;

		// The solution coefficients
		DenseNodeData< Real , Sigs > solution;

		// The average value at the sample positions
		Real isoValue;

		// The density estimator computed from the samples
		DensityEstimator *density;

		// A method that writes the extracted mesh to the streams
		void extractLevelSet( OutputVertexStream< Real , Dim > &vertexStream , OutputFaceStream< Dim-1 > &faceStream , LevelSetExtractionParameters params ) const
		{
			typedef unsigned char AuxData;
			_ExtractLevelSet< false , Real , Dim , FEMSig , AuxData , OutputVertexStream< Real , Dim > , Implicit< Real , Dim , FEMSig > >( IsotropicUIntPack< Dim , FEMSig >() , *this , vertexStream , faceStream , params );
		}
	};

	// Specialized solution information with auxiliary data
	template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
	struct Implicit< Real , Dim , FEMSig , AuxData > : public Implicit< Real , Dim , FEMSig >
	{
		typedef IsotropicUIntPack< Dim , FEMSig > Sigs;

		// The signature of the finite-element used for data extrapolation
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		// The constructor
		Implicit( AuxData zeroAuxData ) : auxData(NULL) , zeroAuxData(zeroAuxData) {}

		// The desctructor
		~Implicit( void ){ delete auxData ; auxData = NULL; }

		// The auxiliary information stored with the oriented vertices
		SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > > *auxData;

		// An instance of "zero" AuxData
		AuxData zeroAuxData;

		// a method for rescaling the contents of the auxiliary data to give interpolation-preference to finer levels
		void weightAuxDataByDepth( Real perLevelScaleFactor )
		{
			auto nodeFunctor = [&]( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > *n )
			{
				ProjectiveData< AuxData , Real >* clr = (*auxData)( n );
				if( clr ) (*clr) *= (Real)pow( (Real)perLevelScaleFactor , Implicit< Real , Dim , FEMSig>::tree.depth( n ) );
			};
			Implicit< Real , Dim , FEMSig>::tree.tree().processNodes( nodeFunctor );
		}

		// A method for writing the extracted mesh to the streams
		void extractLevelSet( OutputVertexWithDataStream< Real , Dim , AuxData > &vertexStream , OutputFaceStream< Dim-1 > &faceStream , LevelSetExtractionParameters params ) const
		{
			_ExtractLevelSet< true , Real , Dim , FEMSig , AuxData , OutputVertexWithDataStream< Real , Dim , AuxData > , Implicit< Real , Dim , FEMSig , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , *this , vertexStream , faceStream , params );
		}
	};

	struct Poisson
	{
		static const unsigned int NormalDegree = 2;							// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
		static const unsigned int DefaultFEMDegree = 1;						// The default finite-element degree (has to be at least 1)
		static const BoundaryType DefaultFEMBoundary = BOUNDARY_NEUMANN;	// The default finite-element boundary type {BOUNDARY_FREE, BOUNDARY_DIRICHLET, BOUNDARY_NEUMANN}
		inline static const float WeightMultiplier = 2.f;					// The default degree-to-point-weight scaling

		template< unsigned int Dim , typename Real >
		struct ConstraintDual
		{
			Real target , weight;
			ConstraintDual( void ) : target(0) , weight(0) {}
			ConstraintDual( Real t , Real w ) : target(t) , weight(w) {}
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p ) const { return CumulativeDerivativeValues< Real , Dim , 0 >( target*weight ); };
		};

		template< unsigned int Dim , typename Real >
		struct SystemDual
		{
			Real weight;
			SystemDual( void ) : weight(0) {}
			SystemDual( Real w ) : weight(w) {}
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
			CumulativeDerivativeValues< double , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< double , Dim , 0 >& dValues ) const { return dValues * weight; };
		};

		template< unsigned int Dim >
		struct SystemDual< Dim , double >
		{
			typedef double Real;
			Real weight;
			SystemDual( void ) : weight(0) {}
			SystemDual( Real w ) : weight(w) {}
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
		};

		template< unsigned int Dim , typename Real >
		struct ValueInterpolationConstraintDual
		{
			typedef VectorTypeUnion< Real , Real > PointSampleData;
			Real vWeight;
			ValueInterpolationConstraintDual( Real v ) : vWeight(v){ }
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim > &p , const VectorTypeUnion< Real , Real >& data ) const 
			{
				Real value = data.template get<0>();
				CumulativeDerivativeValues< Real , Dim , 0 > cdv;
				cdv[0] = value*vWeight;
				return cdv;
			}
		};

		template< unsigned int Dim , typename Real >
		struct ValueInterpolationSystemDual
		{
			CumulativeDerivativeValues< Real , Dim , 0 > weight;
			ValueInterpolationSystemDual( Real v ){ weight[0] = v; }
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Real > &data , const CumulativeDerivativeValues< Real , Dim , 0 > &dValues ) const
			{
				return dValues * weight;
			}
			CumulativeDerivativeValues< double , Dim , 0 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Real > &data , const CumulativeDerivativeValues< double , Dim , 0 > &dValues ) const
			{
				return dValues * weight;
			}
		};

		template< unsigned int Dim >
		struct ValueInterpolationSystemDual< Dim , double >
		{
			typedef double Real;
			Real weight;
			ValueInterpolationSystemDual( void ) : weight(0) {}
			ValueInterpolationSystemDual( Real v ) : weight(v) {}
			CumulativeDerivativeValues< Real , Dim , 0 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Real > &data , const CumulativeDerivativeValues< Real , Dim , 0 > &dValues ) const
			{
				return dValues * weight;
			}
		};

		template< typename Real >
		struct SolutionParameters
		{
			bool verbose;
			bool dirichletErode;
			bool outputDensity;
			bool exactInterpolation;
			bool showResidual;
			Real scale;
			Real confidence;
			Real confidenceBias;
			Real lowDepthCutOff;
			Real width;
			Real pointWeight;
			Real valueInterpolationWeight;
			Real samplesPerNode;
			Real cgSolverAccuracy;
			Real targetValue;
			unsigned int depth;
			unsigned int solveDepth;
			unsigned int baseDepth;
			unsigned int fullDepth;
			unsigned int kernelDepth;
			unsigned int envelopeDepth;
			unsigned int baseVCycles;
			unsigned int iters;
			unsigned int alignDir;

			SolutionParameters( void ) :
				verbose(false) , dirichletErode(false) , outputDensity(false) , exactInterpolation(false) , showResidual(false) ,
				scale((Real)1.1) , confidence((Real)0.) , confidenceBias((Real)0.) , lowDepthCutOff((Real)0.) , width((Real)0.) ,
				pointWeight((Real)0.) , valueInterpolationWeight((Real)0.) , samplesPerNode((Real)1.5) , cgSolverAccuracy((Real)1e-3 ) , targetValue((Real)0.) ,
				depth((unsigned int)8) , solveDepth((unsigned int)-1) , baseDepth((unsigned int)-1) , fullDepth((unsigned int)5) , kernelDepth((unsigned int)-1) ,
				envelopeDepth((unsigned int)-1) , baseVCycles((unsigned int)1) , iters((unsigned int)8) , alignDir(0)
			{}
		};

		template< typename Real , unsigned int Dim >
		struct EnvelopeMesh
		{
			std::vector< Point< Real , Dim > > vertices;
			std::vector< SimplexIndex< Dim-1 , node_index_type > > simplices;
		};

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename ... Other > struct Implicit;

		template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
		static void _Solve( UIntPack< FEMSigs... > , typename std::conditional< HasAuxData , Reconstructor::Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type &implicit , InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh , ValueInterpolationStream< Real , Dim > *valueInterpolationStream );

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename ... Other > struct Implicit;

		template< typename Real , unsigned int Dim , unsigned int FEMSig >
		struct Implicit< Real , Dim , FEMSig > : public Reconstructor::Implicit< Real , Dim , FEMSig >
		{
			Implicit( InputSampleStream< Real , Dim > &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL , ValueInterpolationStream< Real , Dim > *valueInterpolationStream=NULL )
			{
				typedef unsigned char AuxData;
				_Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStream< Real , Dim > >( IsotropicUIntPack< Dim , FEMSig >() , *this , pointStream , params , envelopeMesh , valueInterpolationStream );
			}
		};

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
		struct Implicit< Real , Dim , FEMSig , AuxData > : public Reconstructor::Implicit< Real , Dim , FEMSig , AuxData >
		{
			Implicit( InputSampleWithDataStream< Real , Dim , AuxData > &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL , ValueInterpolationStream< Real , Dim > *valueInterpolationStream=NULL )
				: Reconstructor::Implicit< Real , Dim , FEMSig , AuxData >( pointStream.zero() )
			{
				_Solve< true , Real , Dim , FEMSig , AuxData , InputSampleWithDataStream< Real , Dim , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , *this , pointStream , params , envelopeMesh , valueInterpolationStream );
			}
		};
	};

	struct SSD
	{
		static const unsigned int NormalDegree = 2;									// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
		static const unsigned int DefaultFEMDegree = 2;								// The default finite-element degree (has to be at least 2)
		static const BoundaryType DefaultFEMBoundary = BOUNDARY_NEUMANN;			// The default finite-element boundary type {BOUNDARY_FREE, BOUNDARY_DIRICHLET, BOUNDARY_NEUMANN}
		inline static const double WeightMultipliers[] = { 5e+1f , 5e-4f , 1e-5f };	// The default weights for balancing the value, gradient, and laplacian energy terms

		template< unsigned int Dim , typename ... > struct ConstraintDual;
		template< unsigned int Dim , typename ... > struct SystemDual;

		template< unsigned int Dim , typename Real >
		struct ConstraintDual< Dim , Real >
		{
			Real target , vWeight , gWeight;
			ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const Point< Real , Dim > &n ) const 
			{
				CumulativeDerivativeValues< Real , Dim , 1 > cdv;
				cdv[0] = target*vWeight;
				for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
				return cdv;
			}
		};

		template< unsigned int Dim , typename Real , typename AuxData >
		struct ConstraintDual< Dim , Real , AuxData >
		{
			Real target , vWeight , gWeight;
			ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > &normalAndAuxData ) const 
			{
				Point< Real , Dim > n = normalAndAuxData.template get<0>();
				CumulativeDerivativeValues< Real , Dim , 1 > cdv;
				cdv[0] = target*vWeight;
				for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
				return cdv;
			}
		};

		template< unsigned int Dim , typename Real >
		struct SystemDual< Dim , Real >
		{
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g )
			{
				weight[0] = v;
				for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
			}
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
			CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			};
		};

		template< unsigned int Dim >
		struct SystemDual< Dim , double >
		{
			typedef double Real;
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
		};

		template< unsigned int Dim , typename Real , typename AuxData >
		struct SystemDual< Dim , Real , AuxData >
		{
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g )
			{
				weight[0] = v;
				for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
			}
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
			CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			};
		};

		template< unsigned int Dim , class AuxData >
		struct SystemDual< Dim , double , AuxData >
		{
			typedef double Real;
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
		};

		template< typename Real >
		struct SolutionParameters
		{
			bool verbose;
			bool outputDensity;
			bool exactInterpolation;
			bool showResidual;
			Real scale;
			Real confidence;
			Real confidenceBias;
			Real lowDepthCutOff;
			Real width;
			Real pointWeight;
			Real gradientWeight;
			Real biLapWeight;
			Real samplesPerNode;
			Real cgSolverAccuracy;
			unsigned int depth;
			unsigned int solveDepth;
			unsigned int baseDepth;
			unsigned int fullDepth;
			unsigned int kernelDepth;
			unsigned int baseVCycles;
			unsigned int iters;
			unsigned int alignDir;

			SolutionParameters( void ) :
				verbose(false) , outputDensity(false) , exactInterpolation(false) , showResidual(false) ,
				scale((Real)1.1) , confidence((Real)0.) , confidenceBias((Real)0.) , lowDepthCutOff((Real)0.) , width((Real)0.) ,
				pointWeight((Real)WeightMultipliers[0]) , gradientWeight((Real)WeightMultipliers[1]) , biLapWeight((Real)WeightMultipliers[2]) , samplesPerNode((Real)1.5) , cgSolverAccuracy((Real)1e-3 ) ,
				depth((unsigned int)8) , solveDepth((unsigned int)-1) , baseDepth((unsigned int)-1) , fullDepth((unsigned int)5) , kernelDepth((unsigned int)-1) , alignDir(0) ,
				baseVCycles((unsigned int)1) , iters((unsigned int)8)
			{}

		};

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename ... Other > struct Implicit;

		template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
		static void _Solve( UIntPack< FEMSigs... > , typename std::conditional< HasAuxData , Reconstructor::Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type &implicit , InputSampleStreamType &pointStream , SolutionParameters< Real > params );

		template< typename Real , unsigned int Dim , unsigned int FEMSig >
		struct Implicit< Real , Dim , FEMSig > : public Reconstructor::Implicit< Real , Dim , FEMSig >
		{
			Implicit( InputSampleStream< Real , Dim > &pointStream , SolutionParameters< Real > params )
			{
				typedef unsigned char AuxData;
				_Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStream< Real , Dim > >( IsotropicUIntPack< Dim , FEMSig >() , *this , pointStream , params );
			}
		};

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
		struct Implicit< Real , Dim , FEMSig , AuxData > : public Reconstructor::Implicit< Real , Dim , FEMSig , AuxData >
		{
			Implicit( InputSampleWithDataStream< Real , Dim , AuxData > &pointStream , SolutionParameters< Real > params )
				: Reconstructor::Implicit< Real , Dim , FEMSig , AuxData >( pointStream.zero() )
			{
				_Solve< true , Real , Dim , FEMSig , AuxData , InputSampleWithDataStream< Real , Dim , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , *this , pointStream , params );
			}
		};
	};

	template< class Real , unsigned int Dim , bool ExtendedAxes >
	PointExtent::Extent< Real , Dim , ExtendedAxes > GetExtent( InputDataStream< Point< Real , Dim > > &stream )
	{
		Point< Real , Dim > p;
		PointExtent::Extent< Real , Dim , ExtendedAxes > e;
		while( stream.read( p ) ) e.add( p );
		stream.reset();
		return e;
	}

	template< class Real , unsigned int Dim , bool ExtendedAxes , typename AuxData >
	PointExtent::Extent< Real , Dim , ExtendedAxes > GetExtent( InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , AuxData > > &stream , AuxData d )
	{
		VectorTypeUnion< Real , Point< Real , Dim > , AuxData > p( Point< Real , Dim >() , d );
		PointExtent::Extent< Real , Dim , ExtendedAxes > e;
		while( stream.read( p ) ) e.add( p.template get<0>() );
		stream.reset();
		return e;
	}

	template< class Real , unsigned int Dim , bool ExtendedAxes >
	XForm< Real , Dim+1 > GetPointXForm( InputDataStream< Point< Real , Dim > > &stream , Real scaleFactor , unsigned int dir )
	{
		PointExtent::Extent< Real , Dim , ExtendedAxes > e = GetExtent< Real , Dim , ExtendedAxes >( stream );
		return PointExtent::GetBoundingBoxXForm( e , scaleFactor , dir );
	}

	template< class Real , unsigned int Dim , bool ExtendedAxes , typename AuxData >
	XForm< Real , Dim+1 > GetPointXForm( InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , AuxData > > &stream , AuxData d , Real scaleFactor , unsigned int dir )
	{
		PointExtent::Extent< Real , Dim , ExtendedAxes > e = GetExtent< Real , Dim , ExtendedAxes , AuxData >( stream , d );
		return PointExtent::GetBoundingBoxXForm( e , scaleFactor , dir );
	}

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename OutputVertexStream , typename ImplicitType , unsigned int ... FEMSigs >
	void _ExtractLevelSet( UIntPack< FEMSigs ... > , const ImplicitType &implicit , OutputVertexStream &vertexStream , OutputFaceStream< Dim-1 > &faceStream , LevelSetExtractionParameters params )
	{
		typedef UIntPack< FEMSigs ... > Sigs;
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs ... > >::value , "[ERROR] Signatures don't match" );
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;
		typedef typename ImplicitType::DensityEstimator DensityEstimator;

		if constexpr( Dim==3 )
		{
			Profiler profiler(20);

			std::string statsString;

			if constexpr( HasAuxData )
			{
				typename LevelSetExtractor< Real , Dim , AuxData >::Stats stats;
				TransformedOutputVertexWithDataStream< Real , Dim , AuxData > _vertexStream( implicit.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim , AuxData >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , UIntPack< DataSig >() , implicit.tree , implicit.density , implicit.auxData , implicit.solution , implicit.isoValue , _vertexStream , faceStream , implicit.zeroAuxData , !params.linearFit , params.outputGradients , params.forceManifold , params.polygonMesh , false );
				statsString = stats.toString();
			}
			else
			{
				typename LevelSetExtractor< Real , Dim >::Stats stats;
				TransformedOutputVertexStream< Real , Dim > _vertexStream( implicit.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , implicit.tree , implicit.density , implicit.solution , implicit.isoValue , _vertexStream , faceStream , !params.linearFit , params.outputGradients , params.forceManifold , params.polygonMesh , false );
				statsString = stats.toString();
			}
			if( params.verbose )
			{
				std::cout << "Vertices / Faces: " << vertexStream.size() << " / " << faceStream.size() << std::endl;
				std::cout << statsString << std::endl;
				std::cout << "#            Got Faces: " << profiler << std::endl;
			}
		}
		else if constexpr( Dim==2 )
		{
			Profiler profiler(20);

			std::string statsString;

			if constexpr( HasAuxData )
			{
				typename LevelSetExtractor< Real , Dim , AuxData >::Stats stats;
				TransformedOutputVertexWithDataStream< Real , Dim , AuxData > _vertexStream( implicit.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim , AuxData >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , UIntPack< DataSig >() , implicit.tree , implicit.density , implicit.auxData , implicit.solution , implicit.isoValue , _vertexStream , faceStream , implicit.zeroAuxData , !params.linearFit , params.outputGradients , false );
				statsString = stats.toString();
			}
			else
			{
				typename LevelSetExtractor< Real , Dim >::Stats stats;
				TransformedOutputVertexStream< Real , Dim > _vertexStream( implicit.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , implicit.tree , implicit.density , implicit.solution , implicit.isoValue , _vertexStream , faceStream , !params.linearFit , params.outputGradients , false );
				statsString = stats.toString();
			}
			if( params.verbose )
			{
				std::cout << "Vertices / Faces: " << vertexStream.size() << " / " << faceStream.size() << std::endl;
				std::cout << statsString << std::endl;
				std::cout << "#            Got faces: " << profiler << std::endl;
			}
		}
		else WARN( "Extraction only supported for dimensions 2 and 3" );	}

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
	void Poisson::_Solve( UIntPack< FEMSigs ... > , typename std::conditional< HasAuxData , Reconstructor::Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type &implicit , InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh , ValueInterpolationStream< Real , Dim > *valueInterpolationStream )
	{
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs... > >::value , "[ERROR] Signatures don't match" );
		if( params.valueInterpolationWeight<0 )
		{
			WARN( "Negative value interpolation weight clamped to zero" );
			params.valueInterpolationWeight = 0;
		}
		if( valueInterpolationStream && !params.valueInterpolationWeight ) WARN( "Value interpolation stream provided but interpolation weight is zero" );

		// The signature for the finite-elements representing the auxiliary data (if it's there)
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		///////////////
		// Types --> //
		// The packed finite elements signature
		typedef UIntPack< FEMSigs ... > Sigs;

		// The degrees of the finite elements across the different axes
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;

		// The signature describing the normals elements
		typedef UIntPack< FEMDegreeAndBType< Poisson::NormalDegree , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;

		// Type for tracking sample interpolation
		typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;

		// The finite-element tracking tree node
		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;

		typedef typename FEMTreeInitializer< Dim , Real >::GeometryNodeType GeometryNodeType;

		// The type of the auxiliary information (including the normal)
		typedef typename std::conditional< HasAuxData , VectorTypeUnion< Real , Normal< Real , Dim > , AuxData > , Normal< Real , Dim > >::type NormalAndAuxData;

		// The type describing the sampling density
		typedef typename std::conditional< HasAuxData , Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type::DensityEstimator DensityEstimator;
		// <-- Types //
		///////////////

		NormalAndAuxData zeroNormalAndAuxData;
		if constexpr( HasAuxData ) zeroNormalAndAuxData = NormalAndAuxData( Normal< Real , Dim >() , implicit.zeroAuxData );

		XForm< Real , Dim+1 > modelToUnitCube = XForm< Real , Dim+1 >::Identity();

		Profiler profiler(20);

		size_t pointCount;

		ProjectiveData< Point< Real , 2 > , Real > pointDepthAndWeight;
		std::vector< typename FEMTree< Dim , Real >::PointSample > *valueInterpolationSamples = NULL;
		std::vector< Real > *valueInterpolationSampleData = NULL;
		DenseNodeData< GeometryNodeType , IsotropicUIntPack< Dim , FEMTrivialSignature > > geometryNodeDesignators;
		SparseNodeData< Point< Real , Dim > , NormalSigs > *normalInfo = NULL;
		std::vector< typename FEMTree< Dim , Real >::PointSample > *samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
		std::vector< NormalAndAuxData > *sampleNormalAndAuxData = NULL;

		Real targetValue = params.targetValue;

		// Read in the samples (and auxiliary data)
		{
			profiler.reset();

			pointStream.reset();
			sampleNormalAndAuxData = new std::vector< NormalAndAuxData >();

			modelToUnitCube = params.scale>0 ? GetPointXForm< Real , Dim , true >( pointStream , zeroNormalAndAuxData , params.scale , params.alignDir ) * modelToUnitCube : modelToUnitCube;

			if( params.width>0 )
			{
				// Assuming the transformation is rigid so that the (max) scale can be pulled from the Frobenius norm
				Real maxScale = 0;
				for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) maxScale += modelToUnitCube(i,j) * modelToUnitCube(i,j);
				maxScale = (Real)( 1. / sqrt( maxScale / Dim ) );
				params.depth = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/params.width )/log(2.) ) );
			}
			if( params.solveDepth>params.depth )
			{
				if( params.solveDepth!=-1 ) WARN( "Solution depth cannot exceed system depth: " , params.solveDepth , " <= " , params.depth );
				params.solveDepth = params.depth;
			}
			if( params.fullDepth>params.solveDepth )
			{
				if( params.fullDepth!=-1 ) WARN( "Full depth cannot exceed system depth: " , params.fullDepth , " <= " , params.solveDepth );
				params.fullDepth = params.solveDepth;
			}
			if( params.baseDepth>params.fullDepth )
			{
				if( params.baseDepth!=-1 ) WARN( "Base depth must be smaller than full depth: " , params.baseDepth , " <= " , params.fullDepth );
				params.baseDepth = params.fullDepth;
			}
			if( params.kernelDepth==-1 ) params.kernelDepth = params.depth>2 ? params.depth-2 : 0;
			if( params.kernelDepth>params.depth )
			{
				if( params.kernelDepth!=-1 ) WARN( "Kernel depth cannot exceed system depth: " , params.kernelDepth , " <= " , params.depth );
				params.kernelDepth = params.depth;
			}

			if constexpr( HasAuxData )
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d.template get<0>() /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessData );
			}
			else
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessData );
			}

			implicit.unitCubeToModel = modelToUnitCube.inverse();

			if( params.verbose )
			{
				std::cout << "Input Points / Samples: " << pointCount << " / " << samples->size() << std::endl;
				std::cout << "# Read input into tree: " << profiler << std::endl;
			}
		}

		if( valueInterpolationStream && params.valueInterpolationWeight )
		{
			valueInterpolationSamples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
			valueInterpolationSampleData = new std::vector< Real >();
			// Wrap the point stream in a transforming stream
			TransformedValueInterpolationStream< Real , Dim > _valueInterpolationStream( modelToUnitCube , *valueInterpolationStream );

			// Assign each sample a weight of 1.
			auto ProcessData = []( const Point< Real , Dim > &p , Real &d ){ return (Real)1.; };
			Real zeroValue = (Real)0.;
			typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
			size_t count = FEMTreeInitializer< Dim , Real >::template Initialize< Real >( sid , implicit.tree.spaceRoot() , _valueInterpolationStream , zeroValue , params.depth , *valueInterpolationSamples , *valueInterpolationSampleData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessData );
			if( params.verbose ) std::cout << "Input Interpolation Points / Samples: " << count << " / " << valueInterpolationSamples->size() << std::endl;
		}

		{
			InterpolationInfo *valueInterpolationInfo = NULL;
			DenseNodeData< Real , Sigs > constraints;
			InterpolationInfo *iInfo = NULL;
			int solveDepth = params.depth;

			implicit.tree.resetNodeIndices( 0 , std::make_tuple() );

			// Get the kernel density estimator
			{
				profiler.reset();
				implicit.density = implicit.tree.template setDensityEstimator< 1 , Reconstructor::WeightDegree >( *samples , params.kernelDepth , params.samplesPerNode );
				if( params.verbose ) std::cout << "#   Got kernel density: " << profiler << std::endl;
			}

			// Transform the Hermite samples into a vector field
			{
				profiler.reset();
				normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & ) > ConversionFunction;
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & , Real & ) > ConversionAndBiasFunction;
				if constexpr( HasAuxData )
				{
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				else
				{
					// In this case NormalAndAuxData = Point< Real , Dim >
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				if( params.confidenceBias>0 ) *normalInfo = implicit.tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , implicit.density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionAndBiasFunction );
				else                          *normalInfo = implicit.tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , implicit.density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionFunction );
				ThreadPool::Parallel_for( 0 , normalInfo->size() , [&]( unsigned int , size_t i ){ (*normalInfo)[i] *= (Real)-1.; } );
				if( params.verbose )
				{
					std::cout << "#     Got normal field: " << profiler << std::endl;
					std::cout << "Point depth / Point weight / Estimated measure: " << pointDepthAndWeight.value()[0] << " / " << pointDepthAndWeight.value()[1] << " / " << pointCount*pointDepthAndWeight.value()[1] << std::endl;
				}
			}

			// Get the geometry designators indicating if the space node are interior to, exterior to, or contain the envelope boundary
			if( envelopeMesh )
			{
				profiler.reset();
				{
					// Make the octree complete up to the base depth
					FEMTreeInitializer< Dim , Real >::Initialize( implicit.tree.spaceRoot() , params.baseDepth , []( int , int[] ){ return true; } , implicit.tree.nodeAllocators.size() ?  implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() );

					std::vector< Point< Real , Dim > > vertices( envelopeMesh->vertices.size() );
					for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] = modelToUnitCube * envelopeMesh->vertices[i];
					geometryNodeDesignators = FEMTreeInitializer< Dim , Real >::template GetGeometryNodeDesignators( &implicit.tree.spaceRoot() , vertices , envelopeMesh->simplices , params.baseDepth , params.envelopeDepth , implicit.tree.nodeAllocators , implicit.tree.initializer() );

					// Make nodes in the support of the vector field @{ExactDepth} interior
					if( params.dirichletErode )
					{
						// What to do if we find a node in the support of the vector field
						auto SetScratchFlag = [&]( FEMTreeNode *node )
							{
								if( node )
								{
									while( node->depth()>(int)params.baseDepth ) node = node->parent;
									node->nodeData.setScratchFlag( true );
								}
							};

						std::function< void ( FEMTreeNode * ) > PropagateToLeaves = [&]( const FEMTreeNode *node )
							{
								geometryNodeDesignators[ node ] = GeometryNodeType::INTERIOR;
								if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) PropagateToLeaves( node->children+c );
							};

						// Flags indicating if a node contains a non-zero vector field coefficient
						std::vector< bool > isVectorFieldElement( implicit.tree.nodeCount() , false );

						// Get the set of base nodes
						std::vector< FEMTreeNode * > baseNodes;
						auto nodeFunctor = [&]( FEMTreeNode *node )
							{
								if( node->depth()==params.baseDepth ) baseNodes.push_back( node );
								return node->depth()<(int)params.baseDepth;
							};
						implicit.tree.spaceRoot().processNodes( nodeFunctor );

						std::vector< node_index_type > vectorFieldElementCounts( baseNodes.size() );
						for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCounts[i] = 0;

						// In parallel, iterate over the base nodes and mark the nodes containing non-zero vector field coefficients
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
							{
								auto nodeFunctor = [&]( FEMTreeNode *node )
									{
										Point< Real , Dim > *n = (*normalInfo)( node );
										if( n && Point< Real , Dim >::SquareNorm( *n ) ) isVectorFieldElement[ node->nodeData.nodeIndex ] = true , vectorFieldElementCounts[i]++;
									};
								baseNodes[i]->processNodes( nodeFunctor );
							} );
						size_t vectorFieldElementCount = 0;
						for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCount += vectorFieldElementCounts[i];

						// Get the subset of nodes containing non-zero vector field coefficients and disable the "scratch" flag
						std::vector< FEMTreeNode * > vectorFieldElements;
						vectorFieldElements.reserve( vectorFieldElementCount );
						{
							std::vector< std::vector< FEMTreeNode * > > _vectorFieldElements( baseNodes.size() );
							for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) _vectorFieldElements[i].reserve( vectorFieldElementCounts[i] );
							ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
								{
									auto nodeFunctor = [&]( FEMTreeNode *node )
										{
											if( isVectorFieldElement[ node->nodeData.nodeIndex ] ) _vectorFieldElements[i].push_back( node );
											node->nodeData.setScratchFlag( false );
										};
									baseNodes[i]->processNodes( nodeFunctor );
								} );
							for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) vectorFieldElements.insert( vectorFieldElements.end() , _vectorFieldElements[i].begin() , _vectorFieldElements[i].end() );
						}

						// Set the scratch flag for the base nodes on which the vector field is supported
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] In principal, we should unlock finite elements whose support overlaps the vector field" )
#endif // SHOW_WARNINGS
						implicit.tree.template processNeighboringLeaves< -BSplineSupportSizes< Poisson::NormalDegree >::SupportStart , BSplineSupportSizes< Poisson::NormalDegree >::SupportEnd >( &vectorFieldElements[0] , vectorFieldElements.size() , SetScratchFlag , false );

						// Set sub-trees rooted at interior nodes @ ExactDepth to interior
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ if( baseNodes[i]->nodeData.getScratchFlag() ) PropagateToLeaves( baseNodes[i] ); } );

						// Adjust the coarser node designators in case exterior nodes have become boundary.
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( baseNodes[i] , geometryNodeDesignators ); } );
						FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( &implicit.tree.spaceRoot() , geometryNodeDesignators , params.baseDepth );
					}
				}
				if( params.verbose ) std::cout << "#               Initialized envelope constraints: " << profiler << std::endl;
			}

			if( !params.outputDensity ){ delete implicit.density ; implicit.density = NULL; }
			if constexpr( HasAuxData ) implicit.auxData = new SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > >( implicit.tree.template setExtrapolatedDataField< DataSig , false , Reconstructor::WeightDegree , AuxData >( samples->size() , [&]( size_t i ) -> const typename FEMTree< Dim , Real >::PointSample & { return (*samples)[i]; } , [&]( size_t i ) -> const AuxData & { return (*sampleNormalAndAuxData)[i].template get<1>(); } , (DensityEstimator*)NULL ) );
			delete sampleNormalAndAuxData;

			// Add the interpolation constraints
			if( params.pointWeight>0 )
			{
				profiler.reset();
				if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 > ( implicit.tree , *samples , Poisson::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] ) , Poisson::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] ) , true , false );
				else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( implicit.tree , *samples , Poisson::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] ) , Poisson::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
				if( params.verbose ) std::cout <<  "#Initialized point interpolation constraints: " << profiler << std::endl;
			}

			// Trim the tree and prepare for multigrid
			{
				profiler.reset();
				constexpr int MaxDegree = Poisson::NormalDegree > Degrees::Max() ? Poisson::NormalDegree : Degrees::Max();
				typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
				auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
				auto addNodeFunctor = [&]( int d , const int off[Dim] ){ return d<=(int)params.fullDepth; };
				if constexpr( HasAuxData )
				{
					if( geometryNodeDesignators.size() ) implicit.tree.template finalizeForMultigridWithDirichlet< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(node_index_type)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density , implicit.auxData , &geometryNodeDesignators ) );
					else                                 implicit.tree.template finalizeForMultigrid             < MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor ,                                                                                                                                                                                   std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density , implicit.auxData ) );
				}
				else
				{
					if( geometryNodeDesignators.size() ) implicit.tree.template finalizeForMultigridWithDirichlet< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(node_index_type)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density , &geometryNodeDesignators ) );
					else                                 implicit.tree.template finalizeForMultigrid             < MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor ,                                                                                                                                                                                   std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density ) );
				}

				if( params.verbose ) std::cout << "#       Finalized tree: " << profiler << std::endl;
			}

			// Add the FEM constraints
			{
				profiler.reset();
				constraints = implicit.tree.initDenseNodeData( Sigs() );

				// Add Poisson constraints
				{
					typename FEMIntegrator::template Constraint< Sigs , IsotropicUIntPack< Dim , 1 > , NormalSigs , IsotropicUIntPack< Dim , 0 > , Dim > F;
					unsigned int derivatives2[Dim];
					for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
					typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
					typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
					for( int d=0 ; d<Dim ; d++ )
					{
						unsigned int derivatives1[Dim];
						for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ? 1 : 0;
						F.weights[d][ TensorDerivatives< Derivatives1 >::Index( derivatives1 ) ][ TensorDerivatives< Derivatives2 >::Index( derivatives2 ) ] = 1;
					}
					implicit.tree.addFEMConstraints( F , *normalInfo , constraints , solveDepth );
				}
				if( params.verbose ) std::cout << "#  Set FEM constraints: " << profiler << std::endl;
			}

			// Free up the normal info
			delete normalInfo , normalInfo = NULL;

			if( params.pointWeight>0 )
			{
				profiler.reset();
				implicit.tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "#Set point constraints: " << profiler << std::endl;
			}

			if( valueInterpolationSamples && params.valueInterpolationWeight )
			{
				profiler.reset();
				if( params.exactInterpolation ) valueInterpolationInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , Real , 0 > ( implicit.tree , *valueInterpolationSamples , GetPointer( *valueInterpolationSampleData ) , Poisson::ValueInterpolationConstraintDual< Dim , Real >( params.valueInterpolationWeight ) , Poisson::ValueInterpolationSystemDual< Dim , Real >( params.valueInterpolationWeight ) , true , false );
				else                            valueInterpolationInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , Real , 0 > ( implicit.tree , *valueInterpolationSamples , GetPointer( *valueInterpolationSampleData ) , Poisson::ValueInterpolationConstraintDual< Dim , Real >( params.valueInterpolationWeight ) , Poisson::ValueInterpolationSystemDual< Dim , Real >( params.valueInterpolationWeight ) , true , params.depth , 1 );
				delete valueInterpolationSamples , valueInterpolationSamples = NULL;
				delete valueInterpolationSampleData , valueInterpolationSampleData = NULL;

				implicit.tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( valueInterpolationInfo ) );
				if( params.verbose ) std::cout << "#Set value interpolation constraints: " << profiler << std::endl;
			}

			if( params.verbose ) std::cout << "All Nodes / Active Nodes / Ghost Nodes / Dirichlet Supported Nodes: " << implicit.tree.allNodes() << " / " << implicit.tree.activeNodes() << " / " << implicit.tree.ghostNodes() << " / " << implicit.tree.dirichletElements() << std::endl;
			if( params.verbose ) std::cout << "Memory Usage: " << float( MemoryInfo::Usage())/(1<<20) << " MB" << std::endl;

			// Solve the linear system
			{
				profiler.reset();
				typename FEMTree< Dim , Real >::SolverInfo _sInfo;
				_sInfo.cgDepth = 0 , _sInfo.cascadic = true , _sInfo.vCycles = 1 , _sInfo.iters = params.iters , _sInfo.cgAccuracy = params.cgSolverAccuracy , _sInfo.verbose = params.verbose , _sInfo.showResidual = params.showResidual , _sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , _sInfo.sliceBlockSize = 1;
				_sInfo.baseVCycles = params.baseVCycles;
				typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
				if( valueInterpolationInfo ) implicit.solution = implicit.tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo , valueInterpolationInfo ) );
				else                         implicit.solution = implicit.tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "# Linear system solved: " << profiler << std::endl;
				if( iInfo ) delete iInfo , iInfo = NULL;
				if( valueInterpolationInfo ) delete valueInterpolationInfo , valueInterpolationInfo = NULL;
			}
		}

		// Get the iso-value
		{
			profiler.reset();
			double valueSum = 0 , weightSum = 0;
			typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &implicit.tree , implicit.solution );
			std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
			ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
				{
					ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
					Real w = sample.weight;
					if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
				} );
			for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
			implicit.isoValue = (Real)( valueSum / weightSum );
			if( params.verbose )
			{
				std::cout << "Got average: " << profiler << std::endl;
				std::cout << "Iso-Value: " << implicit.isoValue << " = " << valueSum << " / " << weightSum << std::endl;
			}
		}
		delete samples;
	}

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
	void SSD::_Solve( UIntPack< FEMSigs ... > , typename std::conditional< HasAuxData , Reconstructor::Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type &implicit , InputSampleStreamType &pointStream , SolutionParameters< Real > params )
	{
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs... > >::value , "[ERROR] Signatures don't match" );

		// The signature for the finite-elements representing the auxiliary data (if it's there)
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		///////////////
		// Types --> //
		// The packed finite elements signature
		typedef UIntPack< FEMSigs ... > Sigs;

		// The degrees of the finite elements across the different axes
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;

		// The signature describing the normals elements
		typedef UIntPack< FEMDegreeAndBType< SSD::NormalDegree , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;

		// Type for tracking sample interpolation
		typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 1 > InterpolationInfo;

		// The finite-element tracking tree node
		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;

		// The type of the auxiliary information (including the normal)
		typedef typename std::conditional< HasAuxData , VectorTypeUnion< Real , Normal< Real , Dim > , AuxData > , Normal< Real , Dim > >::type NormalAndAuxData;

		// The type describing the sampling density
		typedef typename std::conditional< HasAuxData , Implicit< Real , Dim , FEMSig , AuxData > , Implicit< Real , Dim , FEMSig > >::type::DensityEstimator DensityEstimator;
		// <-- Types //
		///////////////

		NormalAndAuxData zeroNormalAndAuxData;
		if constexpr( HasAuxData ) zeroNormalAndAuxData = NormalAndAuxData( Normal< Real , Dim >() , implicit.zeroAuxData );

		XForm< Real , Dim+1 > modelToUnitCube = XForm< Real , Dim+1 >::Identity();

		Profiler profiler(20);

		size_t pointCount;

		ProjectiveData< Point< Real , 2 > , Real > pointDepthAndWeight;
		SparseNodeData< Point< Real , Dim > , NormalSigs > *normalInfo = NULL;
		std::vector< typename FEMTree< Dim , Real >::PointSample > *samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
		std::vector< NormalAndAuxData > *sampleNormalAndAuxData = NULL;

		Real targetValue = (Real)0.0;

		// Read in the samples (and auxiliary data)
		{
			profiler.reset();

			pointStream.reset();
			sampleNormalAndAuxData = new std::vector< NormalAndAuxData >();

			modelToUnitCube = params.scale>0 ? GetPointXForm< Real , Dim , true >( pointStream , zeroNormalAndAuxData , params.scale , params.alignDir ) * modelToUnitCube : modelToUnitCube;

			if( params.width>0 )
			{
				// Assuming the transformation is rigid so that the (max) scale can be pulled from the Frobenius norm
				Real maxScale = 0;
				for( unsigned int i=0 ; i<Dim ; i++ ) for( unsigned int j=0 ; j<Dim ; j++ ) maxScale += modelToUnitCube(i,j) * modelToUnitCube(i,j);
				maxScale = (Real)( 1. / sqrt( maxScale / Dim ) );
				params.depth = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/params.width )/log(2.) ) );
			}
			if( params.solveDepth>params.depth )
			{
				if( params.solveDepth!=-1 ) WARN( "Solution depth cannot exceed system depth: " , params.solveDepth , " <= " , params.depth );
				params.solveDepth = params.depth;
			}
			if( params.fullDepth>params.solveDepth )
			{
				if( params.fullDepth!=-1 ) WARN( "Full depth cannot exceed system depth: " , params.fullDepth , " <= " , params.solveDepth );
				params.fullDepth = params.solveDepth;
			}
			if( params.baseDepth>params.fullDepth )
			{
				if( params.baseDepth!=-1 ) WARN( "Base depth must be smaller than full depth: " , params.baseDepth , " <= " , params.fullDepth );
				params.baseDepth = params.fullDepth;
			}
			if( params.kernelDepth==-1 ) params.kernelDepth = params.depth>2 ? params.depth-2 : 0;
			if( params.kernelDepth>params.depth )
			{
				if( params.kernelDepth!=-1 ) WARN( "Kernel depth cannot exceed system depth: " , params.kernelDepth , " <= " , params.depth );
				params.kernelDepth = params.depth;
			}

			if constexpr( HasAuxData )
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d.template get<0>() /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessData );
			}
			else
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , implicit.tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , implicit.tree.nodeAllocators.size() ? implicit.tree.nodeAllocators[0] : NULL , implicit.tree.initializer() , ProcessData );
			}

			implicit.unitCubeToModel = modelToUnitCube.inverse();

			if( params.verbose )
			{
				std::cout << "Input Points / Samples: " << pointCount << " / " << samples->size() << std::endl;
				std::cout << "# Read input into tree: " << profiler << std::endl;
			}
		}
		{
			DenseNodeData< Real , Sigs > constraints;
			InterpolationInfo *iInfo = NULL;
			int solveDepth = params.depth;

			implicit.tree.resetNodeIndices( 0 , std::make_tuple() );

			// Get the kernel density estimator
			{
				profiler.reset();
				implicit.density = implicit.tree.template setDensityEstimator< 1 , Reconstructor::WeightDegree >( *samples , params.kernelDepth , params.samplesPerNode );
				if( params.verbose ) std::cout << "#   Got kernel density: " << profiler << std::endl;
			}

			// Transform the Hermite samples into a vector field
			{
				profiler.reset();
				normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & ) > ConversionFunction;
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & , Real & ) > ConversionAndBiasFunction;
				if constexpr( HasAuxData )
				{
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				else
				{
					// In this case NormalAndAuxData = Point< Real , Dim >
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				if( params.confidenceBias>0 ) *normalInfo = implicit.tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , implicit.density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionAndBiasFunction );
				else                          *normalInfo = implicit.tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , implicit.density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionFunction );
				if( params.verbose )
				{
					std::cout << "#     Got normal field: " << profiler << std::endl;
					std::cout << "Point depth / Point weight / Estimated measure: " << pointDepthAndWeight.value()[0] << " / " << pointDepthAndWeight.value()[1] << " / " << pointCount*pointDepthAndWeight.value()[1] << std::endl;
				}
			}

			if( !params.outputDensity ){ delete implicit.density ; implicit.density = NULL; }
			if constexpr( HasAuxData ) implicit.auxData = new SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > >( implicit.tree.template setExtrapolatedDataField< DataSig , false , Reconstructor::WeightDegree , AuxData >( samples->size() , [&]( size_t i ) -> const typename FEMTree< Dim , Real >::PointSample & { return (*samples)[i]; } , [&]( size_t i ) -> const AuxData & { return (*sampleNormalAndAuxData)[i].template get<1>(); } , (DensityEstimator*)NULL ) );

				// Add the interpolation constraints
				if( params.pointWeight>0 || params.gradientWeight>0 )
				{
					profiler.reset();
						if constexpr( HasAuxData )
						{
							if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , NormalAndAuxData , 1 >( implicit.tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real , AuxData >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real , AuxData >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , false );
							else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , NormalAndAuxData , 1 >( implicit.tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real , AuxData >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real , AuxData >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
						}
						else
						{
							if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , Point< Real , Dim > , 1 >( implicit.tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , false );
							else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , Point< Real , Dim > , 1 >( implicit.tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
						}
					if( params.verbose ) std::cout <<  "#Initialized point interpolation constraints: " << profiler << std::endl;
				}

			delete sampleNormalAndAuxData;


			// Trim the tree and prepare for multigrid
			{
				profiler.reset();
				constexpr int MaxDegree = SSD::NormalDegree > Degrees::Max() ? SSD::NormalDegree : Degrees::Max();
				typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
				auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
				auto addNodeFunctor = [&]( int d , const int off[Dim] ){ return d<=(int)params.fullDepth; };
				if constexpr( HasAuxData ) implicit.tree.template finalizeForMultigrid< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density , implicit.auxData ) );
				else implicit.tree.template finalizeForMultigrid< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , implicit.density ) );

				if( params.verbose ) std::cout << "#       Finalized tree: " << profiler << std::endl;
			}

			// Free up the normal info
			delete normalInfo , normalInfo = NULL;

			if( params.pointWeight>0 || params.gradientWeight>0 )
			{
				profiler.reset();
				constraints = implicit.tree.initDenseNodeData( Sigs() );
				implicit.tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "#Set point constraints: " << profiler << std::endl;
			}

			if( params.verbose ) std::cout << "All Nodes / Active Nodes / Ghost Nodes: " << implicit.tree.allNodes() << " / " << implicit.tree.activeNodes() << " / " << implicit.tree.ghostNodes() << std::endl;
			if( params.verbose ) std::cout << "Memory Usage: " << float( MemoryInfo::Usage())/(1<<20) << " MB" << std::endl;

			// Solve the linear system
			{
				profiler.reset();
				typename FEMTree< Dim , Real >::SolverInfo _sInfo;
				_sInfo.cgDepth = 0 , _sInfo.cascadic = true , _sInfo.vCycles = 1 , _sInfo.iters = params.iters , _sInfo.cgAccuracy = params.cgSolverAccuracy , _sInfo.verbose = params.verbose , _sInfo.showResidual = params.showResidual , _sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , _sInfo.sliceBlockSize = 1;
				_sInfo.baseVCycles = params.baseVCycles;
				typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 2 > > F( { 0. , 0. , (double)params.biLapWeight } );
				implicit.solution = implicit.tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "# Linear system solved: " << profiler << std::endl;
				if( iInfo ) delete iInfo , iInfo = NULL;
			}
		}

		// Get the iso-value
		{
			profiler.reset();
			double valueSum = 0 , weightSum = 0;
			typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &implicit.tree , implicit.solution );
			std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
			ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
				{
					ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
					Real w = sample.weight;
					if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
				} );
			for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
			implicit.isoValue = (Real)( valueSum / weightSum );
			if( params.verbose )
			{
				std::cout << "Got average: " << profiler << std::endl;
				std::cout << "Iso-Value: " << implicit.isoValue << " = " << valueSum << " / " << weightSum << std::endl;
			}
		}
		delete samples;
	}
}


#endif // RECONSTRUCTORS_INCLUDED