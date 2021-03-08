#include "vector3.h"

#include "comm/serialize/serializing_ops.h"

ddl::Writable& operator<< (ddl::Writable &w, const Vector3 &v) {
    return w << v.x << v.y << v.z;
}

ddl::Readable& operator>> (ddl::Readable &r, Vector3 &v) {
    return r >> v.x >> v.y >> v.z;
}
