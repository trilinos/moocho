// ////////////////////////////////////////////////////
// ThyraAccessors.cpp

#include "ThyraAccessors.hpp"
#include "VectorMutableThyra.hpp"
#include "Teuchos_dyn_cast.hpp"

void AbstractLinAlgPack::get_thyra_vector(
	const VectorSpaceThyra                                         &thyra_vec_spc
	,const Vector                                                  &vec
	,Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >    *thyra_vec
	)
{
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
	const VectorMutableThyra *vmtsfc_vec = dynamic_cast<const VectorMutableThyra*>(&vec);
	if(vmtsfc_vec) {
		// We can just grap the const smart pointer to the underlying object 
		*thyra_vec = vmtsfc_vec->thyra_vec();
	}
	else if(vec.space().is_in_core()) {
		// We need to create a temporary copy and then copy the explicit elements
		Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >
			_thyra_vec = ::Thyra::createMember(thyra_vec_spc.thyra_vec_spc());
		// Get explicit views of the elements
		RTOpPack::SubVector vec_sv;
		RTOpPack::MutableSubVector _thyra_vec_sv;
		vec.get_sub_vector( Range1D(), &vec_sv );
		_thyra_vec->getSubVector( Range1D(), &_thyra_vec_sv );
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vec_sv.subDim() != _thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
		// Copy the elements
		for( int i = 0; i < vec_sv.subDim(); ++i )
			_thyra_vec_sv[i] = vec_sv[i];
		// Free/commit the explicit views
		vec.free_sub_vector( &vec_sv );
		_thyra_vec->commitSubVector( &_thyra_vec_sv );
		// Set the output smart pointer
		*thyra_vec = _thyra_vec;
	}
	else {
		TEST_FOR_EXCEPTION(
			true, std::logic_error
			,"AbstractLinAlgPack::get_thyra_vector(...): Error, the vector of concrete type \'"
			<< typeid(vec).name() << "\' is not an incore vector."
			);
	}
}

void AbstractLinAlgPack::free_thyra_vector(
	const VectorSpaceThyra                                         &thyra_vec_spc
	,const Vector                                                  &vec
	,Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >    *thyra_vec
	)
{
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
	*thyra_vec = Teuchos::null;  // This works in both cases above!
}

void AbstractLinAlgPack::get_thyra_vector(
	const VectorSpaceThyra                                         &thyra_vec_spc
	,VectorMutable                                                 *vec
	,Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >          *thyra_vec
	)
{ 
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vec==NULL || thyra_vec==NULL, std::invalid_argument, "Error!" );
#endif
	VectorMutableThyra *vmtsfc_vec = dynamic_cast<VectorMutableThyra*>(vec);
	if(vmtsfc_vec) {
		// We can just directly grap the Thyra vector
		*thyra_vec = vmtsfc_vec->set_uninitialized();
	}
	else if(thyra_vec_spc.is_in_core()) {
		// We need to create a temporary copy and then copy the explicit elements
		Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >
			_thyra_vec = ::Thyra::createMember(thyra_vec_spc.thyra_vec_spc());
		// Get explicit views of the elements
		RTOpPack::SubVector vec_sv;
		vec->get_sub_vector( Range1D(), &vec_sv );
		RTOpPack::MutableSubVector _thyra_vec_sv;
		_thyra_vec->getSubVector( Range1D(), &_thyra_vec_sv );
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vec_sv.subDim() != _thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
		// Copy the elements
		for( int i = 0; i < vec_sv.subDim(); ++i )
			_thyra_vec_sv[i] = vec_sv[i];
		// Free/commit the explicit views
		vec->free_sub_vector( &vec_sv );
		_thyra_vec->commitSubVector( &_thyra_vec_sv );
		// Set the output smart pointer
		*thyra_vec = _thyra_vec;
	}
	else {
		TEST_FOR_EXCEPTION(
			true, std::logic_error
			,"AbstractLinAlgPack::get_thyra_vector(...): Error, the vector of concrete type \'"
			<< typeid(vec).name() << "\' is not an incore vector."
			);
	}
}

void AbstractLinAlgPack::commit_thyra_vector(
	const VectorSpaceThyra                                         &thyra_vec_spc
	,VectorMutable                                                 *vec
	,Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >          *thyra_vec_in
	)
{
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vec==NULL || thyra_vec_in==NULL, std::invalid_argument, "Error!" );
#endif
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >  &thyra_vec = *thyra_vec_in;
	VectorMutableThyra *vmtsfc_vec = dynamic_cast<VectorMutableThyra*>(vec);
	if(vmtsfc_vec) {
		// We can just directly reset the Thyra vector
		vmtsfc_vec->initialize(thyra_vec);
	}
	else if(thyra_vec_spc.is_in_core()) {
		// We need to copy back the temporary
		// Get explicit views of the elements
		RTOpPack::SubVector thyra_vec_sv;
		thyra_vec->getSubVector( Range1D(), &thyra_vec_sv );
		RTOpPack::MutableSubVector vec_sv;
		vec->get_sub_vector( Range1D(), &vec_sv );
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( vec_sv.subDim() != thyra_vec_sv.subDim(), std::logic_error, "Error!" );
#endif
		// Copy the elements
		for( int i = 0; i < vec_sv.subDim(); ++i )
			vec_sv[i] = thyra_vec_sv[i];
		// Free/commit the explicit views
		thyra_vec->freeSubVector( &thyra_vec_sv );
		vec->commit_sub_vector( &vec_sv );
		// Set the output smart pointer
		thyra_vec = Teuchos::null;
	}
	else {
		TEST_FOR_EXCEPTION(	true, std::logic_error, "Should never get here?."	);
	}
}