
#include "base/HarpyInit.h"
#include <map>
#include <vector>
#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <memory>

#include "util/OutputOperators.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

namespace mpi = boost::mpi;
using namespace std;

// Forward declaration of our data class
class DataObject;

// Define a base map type
template<typename Key, typename Val>
using BaseMap = map<Key, vector<Val *>>;

// Sample data class for demonstration purposes
class DataObject {
public:
    int id;
    double value;
    string name;
    vector<string> vs;
    libMesh::RealTensor T;
    libMesh::RealVectorValue RVV;

private:
    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive& ar, const unsigned int /*version*/) {
        ar & id;
        ar & value;
        ar & name;
        ar & vs;
        ar & T;
        ar & RVV;
    }

public:
    DataObject() : id(0), value(0.0), name("") {}
    
    DataObject(int i, double v, const string& n) : id(i), value(v), name(n), T(1,2,3,4,5,6), RVV(1,2,i) {}
    
    void print() const {
        cout << "DataObject[id=" << id << ", value=" << value 
                  << ", name=" << name << "]" ;
        for ( auto & s : vs ) cout << s << "/";
        cout << "|| " << Print(T);
        cout << "|| " << Print(RVV);
        cout << Print(libMesh::Point(1,2,3));
        cout << endl;
    }
};

/**
 *
 *
 *
 *
 */

// Our specialized MPI-enabled map class
template<typename Key, typename Val>
class MPIMap : public BaseMap<Key, Val> {
private:
    mpi::communicator world;
    int rank, size;
    friend class boost::serialization::access;

    // On rank != 0
    void _send() {  world.send(0, 0, *this);  return; }

    // On rank 0
    void _receive( MPIMap<Key,Val> & globalMap ) 
    {
      globalMap = *this;  // add my own data

      for (int i = 1; i < size; ++i) 
      {
        MPIMap<Key,Val> remoteMap;
        world.recv(i, 0, remoteMap);

        for (const auto& [ key, vec ] : remoteMap) 
        {
          auto & tvec = globalMap[key];  // create if non existing
          tvec.insert(tvec.end(), vec.begin(), vec.end());
        }
      }
    }
    
    // For boost
    template<class Archive>
    void serialize(Archive& ar, const unsigned int /*version*/) {
        ar & boost::serialization::base_object<BaseMap<Key,Val>>(*this);
    }

public:
    MPIMap() : BaseMap<Key,Val>(), world(), rank(world.rank()), size(world.size()) {}

    // Main function that collects data from all processes and localizes to rank 0
    void localize_to_one(MPIMap<Key,Val>& globalMap) {
      if (size == 1) { if (rank == 0) globalMap = *this; return; }  // single proc.

      globalMap.clear();
      if (rank != 0) _send(); else _receive( globalMap );
      world.barrier();   // sync
    }
    
};





/**
 *
 *
 *
 *
 *
 *
 */

// Example usage
int main(int argc, char* argv[]) {
    HarpyInit init( argc, argv );
    mpi::communicator world(MPI_COMM_WORLD, mpi::comm_attach);
    
    int rank = world.rank();
    
    // Create local map with some data specific to this process
    MPIMap<string, DataObject> localMap;
    
    // Add some process-specific data
    for (int i = 0; i < 3; i++) {
        int id = rank * 100 + i;
        double value = rank + i / 10.0;
        string name = "Obj_" + to_string(rank) + "_" + to_string(i);
        
        string key = "key_" + to_string(i % 2);  // Create some overlap in keys
        
        DataObject * dd = new DataObject(id, value, name);
        dd->vs.push_back(to_string(rank) + "--");
        dd->vs.push_back(to_string(rank) + "---");
        dd->vs.push_back(to_string(rank) + "----");
        dd->vs.push_back(to_string(rank) + "-----");
        localMap[key].push_back(dd);
    }
    
    // Print local data
    cout << "Process " << rank << " local data:" << endl;
    for (const auto& entry : localMap) {
        cout << "Key: " << entry.first << endl;
        for (const auto& ptr : entry.second) {
            cout << "  ";
            ptr->print();
        }
    }
    
    // Create global map that will contain the combined data on process 0
    MPIMap<string, DataObject> globalMap;
    localMap.localize_to_one(globalMap);
    
    // Print the global data (only process 0 will have it)
    if (rank == 0) {
        cout << "\nGlobal data on process 0:" << endl;
        for (const auto& entry : globalMap) {
            cout << "Key: " << entry.first << endl;
            for (const auto& ptr : entry.second) {
                cout << "  ";
                ptr->print();
            }
        }
    }
    
    return 0;
}



namespace boost { namespace serialization {
/** **/
template<class Archive, typename T>
void serialize(Archive & ar, libMesh::TensorValue<T> & tensor, const unsigned int /*version*/)
{
  for (unsigned int i = 0; i < 3; ++i) 
  for (unsigned int j = 0; j < 3; ++j) 
    ar & tensor(i,j);
} 
/** **/
template<class Archive, typename T>
void serialize(Archive & ar, libMesh::VectorValue<T> & vector, const unsigned int /*version*/)
{ for (uint i=0; i<3; i++ ) ar & vector(i); }

} } // namespace boost
