// ///////////////////////////////////////////////////////
// Serializable.hpp

#ifndef SERIALIZATIONPACK_SERIALIZABLE_HPP
#define SERIALIZATIONPACK_SERIALIZABLE_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace SerializationPack {

///
/** Mixin interface for objects that can be serialized to and from a stream.
 *
 * Todo: Finish documentation!
 */
class Serializable {
public:
	
	///
	virtual ~Serializable() {}

	///
	/** Serialize the object to a stream.
	 *
	 * Todo: Finish documentation!
	 */
	virtual void serialize( std::ostream &out ) const = 0;

	///
	/** Un-serialize the object from a stream.
	 *
	 * Todo: Finish documentation!
	 */
	virtual void unserialize( std::istream &in ) = 0;

};

} // namespace SerializationPack

#endif // SERIALIZATIONPACK_SERIALIZABLE_HPP