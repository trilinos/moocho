// //////////////////////////////////////////////////////////////////////
// SpVectorOp.cpp

#include "SparseLinAlgPack/include/SpVectorOp.h"

namespace {
// Setup some template classes to check at complile time that
// the layout of SpVectorSlice::element_type is proper.
template<int N, class T>
class assert_compile_time {
public:
assert_compile_time()
{
	// This should not compile if instantiated with a type T that
	// is not an integer.  However, if the compiler checks this
	// function without instantiating it, it can not cause an error
	// because it does not know the type of T to see if the
	// conversion is legal or not.
	T d;
	static_cast<int*>(d);
}
};
// Template specialization for no error
template<>
class assert_compile_time<0,double> {
public:
assert_compile_time()
{}
};
// Validate that there is an integer stride between values
assert_compile_time<
    ((int)sizeof(SparseLinAlgPack::SpVectorSlice::element_type)
	 % (int)sizeof(LinAlgPack::value_type))
	, double
	>
    validate_value_stride;
// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(SparseLinAlgPack::SpVectorSlice::element_type)
	 % (int)sizeof(LinAlgPack::indice_type))
	, double
	>
    validate_index_stride;
} // end namespace

void SparseLinAlgPack::add_elements( SpVector* sv_lhs, const VectorSlice& vs_rhs, size_type offset )
{
	typedef SpVector::element_type ele_t;
	const bool assume_sorted = !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() );
	VectorSlice::const_iterator
		itr = vs_rhs.begin();
	for( size_type i = 1; i <= vs_rhs.size(); ++i )
		sv_lhs->add_element( ele_t( i + offset, *itr++ ) );
	sv_lhs->assume_sorted(assume_sorted);
}

void SparseLinAlgPack::add_elements( SpVector* sv_lhs, const SpVectorSlice& sv_rhs, size_type offset )
{
	typedef SpVector::element_type ele_t;
	const bool assume_sorted = ( !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() ) )
		&& ( !sv_rhs.nz() || ( sv_rhs.nz() || sv_rhs.is_sorted() ) );
	for( SpVectorSlice::const_iterator itr = sv_rhs.begin(); itr != sv_rhs.end(); ++itr )
		sv_lhs->add_element( ele_t( itr->indice() + sv_rhs.offset() + offset, itr->value() ) );
	sv_lhs->assume_sorted(assume_sorted);
}