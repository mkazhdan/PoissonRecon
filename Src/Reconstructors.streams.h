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

#ifndef RECONSTRUCTORS_STREAMS_INCLUDED
#define RECONSTRUCTORS_STREAMS_INCLUDED

#include <mutex>

// Basic types
template< typename Real , unsigned int Dim >
using BaseSample = VectorTypeUnion< Real , Position< Real , Dim > , Normal< Real , Dim > >;

template< typename Real , unsigned int Dim , typename Data >
using BaseSampleWithData = VectorTypeUnion< Real , Position< Real , Dim > , VectorTypeUnion< Real , Normal< Real , Dim > , Data > >;

template< typename Real , unsigned int Dim >
using BaseVertex = VectorTypeUnion< Real , Position< Real , Dim > , Normal< Real , Dim > , Real >;

template< typename Real , unsigned int Dim , typename Data >
using BaseVertexWithData = VectorTypeUnion< Real , Position< Real , Dim > , Gradient< Real , Dim > , Real , Data >;

using Edge = std::pair< node_index_type , node_index_type >;
using Polygon = std::vector< node_index_type >;

template< unsigned int FaceDim > using Face = std::conditional_t< FaceDim==2 , Polygon , std::conditional_t< FaceDim==1 , Edge , void * > >;

// Basic streams
template< typename Real , unsigned int Dim                 > using BaseInputSampleStream         = InputDataStream< BaseSample        < Real , Dim >        >;
template< typename Real , unsigned int Dim , typename Data > using BaseInputSampleWithDataStream = InputDataStream< BaseSampleWithData< Real , Dim , Data > >;

template< typename Real , unsigned int Dim                 > using BaseOutputVertexStream         = OutputDataStream< BaseVertex        < Real , Dim >        >;
template< typename Real , unsigned int Dim , typename Data > using BaseOutputVertexWithDataStream = OutputDataStream< BaseVertexWithData< Real , Dim , Data > >;

template< unsigned int FaceDim > using  InputFaceStream =  InputDataStream< Face< FaceDim > >;
template< unsigned int FaceDim > using OutputFaceStream = OutputDataStream< Face< FaceDim > >;

/////////////////////////////////////
// Value Interpolation Data Stream //
/////////////////////////////////////
template< typename Real , unsigned int Dim >
struct ValueInterpolationStream : public InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , Real > >
{
	// Functionality to reset the stream to the start
	virtual void reset( void ) = 0;

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	virtual bool base_read( Position< Real , Dim > &p , Real &v ) = 0;
	// Implementation of InputDataStream::read
	bool base_read( VectorTypeUnion< Real , Point< Real , Dim > , Real > &s ){ return base_read( s.template get<0>() , s.template get<1>() ); }
};

/////////////////////////////////////////////////
// Transformed Value Interpolation Data Stream //
/////////////////////////////////////////////////
template< typename Real , unsigned int Dim >
struct TransformedValueInterpolationStream : public ValueInterpolationStream< Real , Dim >
{
	// A constructor initialized with the transformation to be applied to the samples, and a sample stream
	TransformedValueInterpolationStream( XForm< Real , Dim+1 > xForm , ValueInterpolationStream< Real , Dim > &stream ) : _stream(stream) , _xForm(xForm) {}

	// Functionality to reset the stream to the start
	void reset( void ){ _stream.reset(); }

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	bool base_read( Position< Real , Dim > &p , Real &v )
	{
		VectorTypeUnion< Real , Point< Real , Dim > , Real > s;
		bool ret = _stream.read( s );
		if( ret ) p = _xForm * s.template get<0>() , v = s.template get<1>();
		return ret;
	}

protected:
	// A reference to the underlying stream
	ValueInterpolationStream< Real , Dim > &_stream;

	// The affine transformation to be applied to the positions
	XForm< Real , Dim+1 > _xForm;
};

///////////////////////////
// Oriented Point Stream //
///////////////////////////
template< typename Real , unsigned int Dim >
struct InputSampleStream : public BaseInputSampleStream< Real , Dim >
{
	// Functionality to reset the stream to the start
	virtual void reset( void ) = 0;

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	virtual bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n ) = 0;
	// Implementation of InputDataStream::read
	bool base_read( BaseSample< Real , Dim > &s ){ return base_read( s.template get<0>() , s.template get<1>() ); }
};

///////////////////////////////////
// Oriented Point w/ Data Stream //
///////////////////////////////////
template< typename Real , unsigned int Dim , typename Data >
struct InputSampleWithDataStream : public BaseInputSampleWithDataStream< Real , Dim , Data >
{
	// A constructor initialized with an instance of "zero" data
	InputSampleWithDataStream( Data zero ) : _zero(zero) {}

	// Functionality to reset the stream to the start
	virtual void reset( void ) = 0;

	// Returns the zero instance
	const Data &zero( void ) const{ return _zero; }

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	virtual bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d ) = 0;
	bool base_read( BaseSampleWithData< Real , Dim , Data > &s ){ return base_read( s.template get<0>() , s.template get<1>().template get<0>() , s.template get<1>().template get<1>() ); }

	// An instance of "zero" data
	Data _zero;
};


///////////////////////////////////////
// Transformed Oriented Point Stream //
///////////////////////////////////////
#ifdef DE_VIRTUALIZE_INPUT
template< typename Real , unsigned int Dim , typename InputStream >
#else // !DE_VIRTUALIZE_INPUT
template< typename Real , unsigned int Dim >
#endif // DE_VIRTUALIZE_INPUT
struct TransformedInputSampleStream : public InputSampleStream< Real , Dim >
{
#ifdef DE_VIRTUALIZE_INPUT
	static_assert( std::is_base_of< InputSampleStream< Real , Dim > , InputStream >::value , "[ERROR] Unexpected stream type" );
#endif // DE_VIRTUALIZE_INPUT
	// A constructor initialized with the transformation to be applied to the samples, and a sample stream
#ifdef DE_VIRTUALIZE_INPUT
	TransformedInputSampleStream( XForm< Real , Dim+1 > xForm , InputStream &stream ) : _stream(stream) , _positionXForm(xForm)
#else // !DE_VIRTUALIZE_INPUT
	TransformedInputSampleStream( XForm< Real , Dim+1 > xForm , InputSampleStream< Real , Dim > &stream ) : _stream(stream) , _positionXForm(xForm)
#endif // DE_VIRTUALIZE_INPUT
	{
		_normalXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( fabs( xForm.determinant() ) , 1./Dim );
	}

	// Functionality to reset the stream to the start
	void reset( void ){ _stream.reset(); }

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n )
	{
		BaseSample< Real , Dim > s;
		bool ret = _stream.read( s );
		if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>();
		return ret;
	}

protected:
	// A reference to the underlying stream
#ifdef DE_VIRTUALIZE_INPUT
	InputStream &_stream;
#else // !DE_VIRTUALIZE_INPUT
	InputSampleStream< Real , Dim > &_stream;
#endif // DE_VIRTUALIZE_INPUT

	// The affine transformation to be applied to the positions
	XForm< Real , Dim+1 > _positionXForm;

	// The linear transformation to be applied to the normals
	XForm< Real , Dim > _normalXForm;
};

///////////////////////////////////////////////
// Transformed Oriented Point w/ Data Stream //
///////////////////////////////////////////////
#ifdef DE_VIRTUALIZE_INPUT
template< typename Real , unsigned int Dim , typename Data , typename InputStream >
#else // !DE_VIRTUALIZE_INPUT
template< typename Real , unsigned int Dim , typename Data >
#endif // DE_VIRTUALIZE_INPUT
struct TransformedInputSampleWithDataStream : public InputSampleWithDataStream< Real , Dim , Data >
{
#ifdef DE_VIRTUALIZE_INPUT
	static_assert( std::is_base_of< InputSampleWithDataStream< Real , Dim , Data > , InputStream >::value , "[ERROR] Unexpected stream type" );
#endif // DE_VIRTUALIZE_INPUT
	// A constructor initialized with an instance of "zero" data
#ifdef DE_VIRTUALIZE_INPUT
	TransformedInputSampleWithDataStream( XForm< Real , Dim+1 > xForm , InputStream &stream ) : InputSampleWithDataStream< Real , Dim , Data >( stream.zero() ) , _stream(stream) , _positionXForm(xForm)
#else // !DE_VIRTUALIZE_INPUT
	TransformedInputSampleWithDataStream( XForm< Real , Dim+1 > xForm , InputSampleWithDataStream< Real , Dim , Data > &stream ) : InputSampleWithDataStream< Real , Dim , Data >( stream.zero() ) , _stream(stream) , _positionXForm(xForm)
#endif // DE_VIRTUALIZE_INPUT
	{
		_normalXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
	}

	// Functionality to reset the stream to the start
	void reset( void ){ _stream.reset(); }

	// Functionality to extract the next position/normal pair.
	// The method returns true if there was another point in the stream to read, and false otherwise
	bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d )
	{
		BaseSampleWithData< Real , Dim , Data > s( Position< Real , Dim >() , VectorTypeUnion< Real , Normal< Real , Dim > , Data  >( Normal< Real , Dim >() , _stream.zero() ) );
		bool ret = _stream.read( s );
		if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>().template get<0>() , d = s.template get<1>().template get<1>();
		return ret;
	}

protected:
	// A reference to the underlying stream
#ifdef DE_VIRTUALIZE_INPUT
	InputStream &_stream;
#else // !DE_VIRTUALIZE_INPUT
	InputSampleWithDataStream< Real , Dim , Data > &_stream;
#endif // DE_VIRTUALIZE_INPUT

	// The affine transformation to be applied to the positions
	XForm< Real , Dim+1 > _positionXForm;

	// The linear transformation to be applied to the normals
	XForm< Real , Dim > _normalXForm;
};


///////////////////
// Vertex Stream //
///////////////////
template< typename Real , unsigned int Dim >
struct OutputVertexStream : public BaseOutputVertexStream< Real , Dim >
{
	// Need to provide access to base write for counter support
	using BaseOutputVertexStream< Real , Dim >::write;

	// Functionality to insert the next vertex
	virtual void base_write( Position< Real , Dim > p , Gradient< Real , Dim > g , Real w ) = 0;
	void base_write( const BaseVertex< Real , Dim > &v ){ base_write( v.template get<0>() , v.template get<1>() , v.template get<2>() ); }
};

///////////////////////////
// Vertex w/ Data Stream //
///////////////////////////
template< typename Real , unsigned int Dim , typename Data >
struct OutputVertexWithDataStream : public BaseOutputVertexWithDataStream< Real , Dim , Data >
{
	// Need to provide access to base write for counter support
	using BaseOutputVertexWithDataStream< Real , Dim , Data >::write;

	// Functionality to insert the next vertex
	virtual void base_write( Position< Real , Dim > p , Gradient< Real , Dim > g , Real w , Data d ) = 0;
	void base_write( const BaseVertexWithData< Real , Dim , Data > &v ){ return base_write( v.template get<0>() , v.template get<1>() , v.template get<2>() , v.template get<3>() ); }
};

///////////////////////////////
// Transformed Vertex Stream //
///////////////////////////////
template< typename Real , unsigned int Dim >
struct TransformedOutputVertexStream : public OutputVertexStream< Real , Dim >
{
	// A constructor initialized with the transformation to be applied to the samples, and a sample stream
	TransformedOutputVertexStream( XForm< Real , Dim+1 > xForm , OutputVertexStream< Real , Dim > &stream ) : _stream(stream) , _positionXForm(xForm)
	{
		_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
	}

	// Need to write the union to ensure that the counter gets set
	void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w ){ _stream.write( BaseVertex< Real , Dim >( _positionXForm * p , _gradientXForm * g , w ) ); }

protected:
	// A reference to the underlying stream
	OutputVertexStream< Real , Dim > &_stream;

	// The affine transformation to be applied to the positions
	XForm< Real , Dim+1 > _positionXForm;

	// The linear transformation to be applied to the normals
	XForm< Real , Dim > _gradientXForm;
};


///////////////////////////////////////
// Transformed Vertex w/ Data Stream //
///////////////////////////////////////
template< typename Real , unsigned int Dim , typename Data >
struct TransformedOutputVertexWithDataStream : public OutputVertexWithDataStream< Real , Dim , Data >
{
	// A constructor initialized with the transformation to be applied to the samples, and a sample stream
	TransformedOutputVertexWithDataStream( XForm< Real , Dim+1 > xForm , OutputVertexWithDataStream< Real , Dim , Data > &stream ) : _stream(stream) , _positionXForm(xForm)
	{
		_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
	}

	void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d ){ _stream.write( BaseVertexWithData< Real , Dim , Data >( _positionXForm * p , _gradientXForm * g , w , d ) ); }

protected:
	// A reference to the underlying stream
	OutputVertexWithDataStream< Real , Dim , Data > &_stream;

	// The affine transformation to be applied to the positions
	XForm< Real , Dim+1 > _positionXForm;

	// The linear transformation to be applied to the normals
	XForm< Real , Dim > _gradientXForm;
};

///////////////////////////////////////////
// A wrapper class to write out vertices //
///////////////////////////////////////////
template< typename Real , unsigned int Dim , typename Vertex >
struct OutputVertexStreamWrapper : public OutputVertexStream< Real , Dim >
{
	virtual void set( Vertex &out , const BaseVertex< Real , Dim > &in ) = 0;

	OutputVertexStreamWrapper( OutputDataStream< Vertex > &stream , Vertex out ) : _stream(stream) , _out(out) {}

	void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
	{
		_in.template get<0>() = p;
		_in.template get<1>() = g;
		_in.template get<2>() = w;
		set( _out , _in );
		_stream.write( _out );
	}
protected:
	OutputDataStream< Vertex > &_stream;
	BaseVertex< Real , Dim > _in;
	Vertex _out;
};

template< typename Real , unsigned int Dim , typename Data , typename Vertex >
struct OutputVertexWithDataStreamWrapper : public OutputVertexWithDataStream< Real , Dim , Data >
{
	virtual void set( Vertex &out , const BaseVertexWithData< Real , Dim , Data > &in ) = 0;

	OutputVertexWithDataStreamWrapper( OutputDataStream< Vertex > &stream , BaseVertexWithData< Real , Dim , Data > in , Vertex out ) : _stream(stream) , _in(in) , _out(out) {}

	void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
	{
		_in.template get<0>() = p;
		_in.template get<1>() = g;
		_in.template get<2>() = w;
		_in.template get<3>() = d;
		set( _out , _in );
		_stream.write( _out );
	}
protected:
	OutputDataStream< Vertex > &_stream;
	BaseVertexWithData< Real , Dim , Data > _in;
	Vertex _out;
};


//////////////////////////////////////////////////////////////////////////////
// Output and the input face stream, backed either by a file or by a vector //
//////////////////////////////////////////////////////////////////////////////
// [WARNING] These assume that the stream starts as write-only and after the reset method is invoked, the stream becomes read-only.

template< unsigned int FaceDim >
struct OutputInputFaceStream : public OutputFaceStream< FaceDim > , public InputFaceStream< FaceDim >
{
	// The streams for communicating the information
	InputFaceStream < FaceDim > * inStream;
	OutputFaceStream< FaceDim > *outStream;

	void reset( void ){ inStream->reset(); }
	bool base_read( Face< FaceDim > &f ){ return inStream->read(f); }
	bool base_read( unsigned int t , Face< FaceDim > &f ){ return inStream->read(t,f); }
	void base_write( const Face< FaceDim > &f ){ outStream->write(f); }
	void base_write( unsigned int t , const Face< FaceDim > &f ){ outStream->write(t,f); }

	OutputInputFaceStream( bool inCore , bool multi )
	{
		size_t sz = std::thread::hardware_concurrency();

		_backingVector = NULL;
		_backingVectors.resize( sz , NULL );

		_backingFile = NULL;
		_backingFiles.resize( sz , NULL );

		_inStreams.resize( sz , NULL );
		_outStreams.resize( sz , NULL );

		if( inCore )
		{
			if( multi )
			{
				for( unsigned int i=0 ; i<sz ; i++ )
				{
					_backingVectors[i] = new std::vector< Face< FaceDim > >();
					_inStreams[i] = new VectorBackedInputDataStream< Face< FaceDim > >( *_backingVectors[i] );
					_outStreams[i] = new VectorBackedOutputDataStream< Face< FaceDim > >( *_backingVectors[i] );
				}
				inStream = new MultiInputDataStream< Face< FaceDim > >( _inStreams );
				outStream = new MultiOutputDataStream< Face< FaceDim > >( _outStreams );
			}
			else
			{
				_backingVector = new std::vector< Face< FaceDim > >();
				inStream = new VectorBackedInputDataStream< Face< FaceDim > >( *_backingVector );
				outStream = new VectorBackedOutputDataStream< Face< FaceDim > >( *_backingVector );
			}
		}
		else
		{
			if( multi )
			{
				for( unsigned int i=0 ; i<sz ; i++ )
				{
					_backingFiles[i] = new FileBackedReadWriteStream::FileDescription( NULL );
					_inStreams[i] = new FileBackedInputDataStream< Face< FaceDim > >( _backingFiles[i]->fp );
					_outStreams[i] = new FileBackedOutputDataStream< Face< FaceDim > >( _backingFiles[i]->fp );
				}
				inStream = new MultiInputDataStream< Face< FaceDim > >( _inStreams );
				outStream = new MultiOutputDataStream< Face< FaceDim > >( _outStreams );
			}
			else
			{
				_backingFile = new FileBackedReadWriteStream::FileDescription( NULL );
				inStream = new FileBackedInputDataStream< Face< FaceDim > >( _backingFile->fp );
				outStream = new FileBackedOutputDataStream< Face< FaceDim > >( _backingFile->fp );
			}
		}
	}

	~OutputInputFaceStream( void )
	{
		size_t sz = std::thread::hardware_concurrency();

		delete _backingVector;
		delete _backingFile;

		for( unsigned int i=0 ; i<sz ; i++ )
		{
			delete _backingVectors[i];
			delete _backingFiles[i];
			delete  _inStreams[i];
			delete _outStreams[i];
		}

		delete  inStream;
		delete outStream;
	}
protected:
	std::vector< Face< FaceDim > > *_backingVector;
	FileBackedReadWriteStream::FileDescription *_backingFile;
	std::vector< std::vector< Face< FaceDim > > * > _backingVectors;
	std::vector< FileBackedReadWriteStream::FileDescription * > _backingFiles;
	std::vector<  InputDataStream< Face< FaceDim > > * >  _inStreams;
	std::vector< OutputDataStream< Face< FaceDim > > * > _outStreams;
};

template< typename Factory >
struct OutputInputFactoryTypeStream : public OutputDataStream< typename Factory::VertexType > , public InputDataStream< typename Factory::VertexType >
{
	typedef typename Factory::VertexType Vertex;
	// The streams for communicating the information
	InputDataStream < Vertex > * inStream;
	OutputDataStream< Vertex > *outStream;

	void reset( void ){ inStream->reset(); }
	void base_write( const Vertex &v ){ outStream->write( v ); }
	bool base_read( Vertex &v ){ return inStream->read( v ); }

	OutputInputFactoryTypeStream( Factory &factory , bool inCore , bool multi )
	{
		size_t sz = std::thread::hardware_concurrency();

		_backingVector = NULL;
		_backingVectors.resize( sz , NULL );

		_backingFile = NULL;
		_backingFiles.resize( sz , NULL );

		_inStreams.resize( sz , NULL );
		_outStreams.resize( sz , NULL );

		if( inCore )
		{
			if( multi )
			{
				for( unsigned int i=0 ; i<sz ; i++ )
				{
					_backingVectors[i] = new std::vector< Vertex >();
					_inStreams[i] = new VectorBackedInputDataStream< Vertex >( *_backingVectors[i] );
					_outStreams[i] = new VectorBackedOutputDataStream< Vertex >( *_backingVectors[i] );
				}
				inStream = new MultiInputDataStream< Vertex >( _inStreams );
				outStream = new MultiOutputDataStream< Vertex >( _outStreams );
			}
			else
			{
				_backingVector = new std::vector< Vertex >();

				inStream = new VectorBackedInputDataStream< Vertex >( *_backingVector );
				outStream = new VectorBackedOutputDataStream< Vertex >( *_backingVector );
			}
		}
		else
		{
			if( multi )
			{
				for( unsigned int i=0 ; i<sz ; i++ )
				{
					_backingFiles[i] = new FileBackedReadWriteStream::FileDescription( NULL );
					_inStreams[i] = new FileBackedInputFactoryTypeStream< Factory >( _backingFiles[i]->fp , factory );
					_outStreams[i] = new FileBackedOutputFactoryTypeStream< Factory >( _backingFiles[i]->fp , factory  );
				}
				inStream = new MultiInputDataStream< Vertex >( _inStreams );
				outStream = new MultiOutputDataStream< Vertex >( _outStreams );
			}
			else
			{
				_backingFile = new FileBackedReadWriteStream::FileDescription( NULL );
				inStream = new FileBackedInputFactoryTypeStream< Factory >( _backingFile->fp , factory );
				outStream = new FileBackedOutputFactoryTypeStream< Factory >( _backingFile->fp , factory );
			}
		}
	}

	~OutputInputFactoryTypeStream( void )
	{
		size_t sz = std::thread::hardware_concurrency();

		delete _backingVector;
		delete _backingFile;

		for( unsigned int i=0 ; i<sz ; i++ )
		{
			delete _backingVectors[i];
			delete _backingFiles[i];
			delete  _inStreams[i];
			delete _outStreams[i];
		}

		delete  inStream;
		delete outStream;
	}
protected:
	std::vector< Vertex > *_backingVector;
	FileBackedReadWriteStream::FileDescription *_backingFile;
	std::vector< std::vector< Vertex > * >_backingVectors;
	std::vector< FileBackedReadWriteStream::FileDescription * > _backingFiles;
	std::vector<  InputDataStream< Vertex > * >  _inStreams;
	std::vector< OutputDataStream< Vertex > * > _outStreams;
};

template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity >
struct OutputVertexInfo
{
	using Factory =
		typename std::conditional
		<
			HasGradients ,
			typename std::conditional
			<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > >
			>::type ,
			typename std::conditional
			<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::PositionFactory< Real , Dim >
			>::type
		>::type;
	using Vertex = typename Factory::VertexType;

	static Factory GetFactory( void ){ return Factory(); }

	struct StreamWrapper : public Reconstructor::OutputVertexStreamWrapper< Real , Dim , Vertex >
	{
		StreamWrapper( OutputDataStream< Vertex > &stream , Vertex out ) :
			Reconstructor::OutputVertexStreamWrapper< Real , Dim , Vertex >( stream , out ){}
		void set( Vertex &out , const Reconstructor::BaseVertex< Real , Dim > &in )
		{
			if constexpr( HasGradients || HasDensity )
			{
				out.template get<0>() = in.template get<0>();
				if constexpr( HasGradients )
				{
					out.template get<1>() = in.template get<1>();
					if constexpr( HasDensity ) out.template get<2>() = in.template get<2>();
				}
				else
				{
					if constexpr( HasDensity ) out.template get<1>() = in.template get<2>();
				}
			}
			else out = in.template get<0>();
		}
	};
};

template< typename Real , unsigned int Dim , typename AuxDataFactory , bool HasGradients , bool HasDensity >
struct OutputVertexWithDataInfo
{
	using Factory =
		typename std::conditional
		<
			HasGradients ,
			typename std::conditional
			<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , AuxDataFactory >
			>::type ,
			typename std::conditional
			<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , AuxDataFactory >
			>::type
		>::type;
	using AuxData = typename AuxDataFactory::VertexType;

	using _Vertex = VectorTypeUnion< Real , Point< Real , Dim > , Point< Real , Dim > , Real , typename AuxDataFactory::VertexType >;
	using Vertex = typename Factory::VertexType;

	static Factory GetFactory( AuxDataFactory auxDataFactory )
	{
		if constexpr( HasGradients )
		{
			if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() , VertexFactory::ValueFactory< Real >() , auxDataFactory );
			else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() ,                                         auxDataFactory );
		}
		else
		{
			if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                VertexFactory::ValueFactory< Real >() , auxDataFactory );
			else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                                                        auxDataFactory );
		}
	}

	struct StreamWrapper : public Reconstructor::OutputVertexWithDataStreamWrapper< Real , Dim , AuxData , Vertex >
	{
		StreamWrapper( OutputDataStream< Vertex > &stream , Vertex out ) :
			Reconstructor::OutputVertexWithDataStreamWrapper< Real , Dim , AuxData , Vertex >( stream , _Vertex() , out ){}

		void set( Vertex &out , const Reconstructor::BaseVertexWithData< Real , Dim , AuxData > &in )
		{
			out.template get<0>() = in.template get<0>();
			if constexpr( HasGradients )
			{
				out.template get<1>() = in.template get<1>();
				if constexpr( HasDensity ) out.template get<2>() = in.template get<2>() , out.template get<3>() = in.template get<3>();
				else out.template get<2>() = in.template get<3>();
			}
			else
			{
				if constexpr( HasDensity ) out.template get<1>() = in.template get<2>() , out.template get<2>() = in.template get<3>();
				else out.template get<1>() = in.template get<3>();
			}
		}
	};
};

#endif // RECONSTRUCTORS_STREAMS_INCLUDED