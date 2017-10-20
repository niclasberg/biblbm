#ifndef COMMUNICATIONBUFFER_H_
#define COMMUNICATIONBUFFER_H_
#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "geometry.h"
#include <map>
#include <vector>

namespace plb {

namespace fsi {

class CommunicationBuffer {
public:
	typedef std::map<plint, std::vector<char> > BufferType;
	typedef BufferType::iterator BufferIterator;

	CommunicationBuffer()
	: transfer_in_progress(false), proc_list(), send_requests(0), recv_requests(0),
	  recv_request_count(0), send_request_count(0)
	{

	}

	CommunicationBuffer(const std::vector<plint> & proc_list)
	: transfer_in_progress(false), proc_list(proc_list), send_requests(0), recv_requests(0)
	{
		allocate_request_arrays();
	}

	virtual ~CommunicationBuffer() { free_request_arrays(); }

	void set_proc_list(const std::vector<plint> & proc_list) { this->proc_list = proc_list; allocate_request_arrays(); }

	void send_and_receive_no_wait();
	void finalize_send_and_receive();

	void clear_send_buffer();

	pluint bytes_received() const { return recv_buffer.size(); }
	char * recv_buffer_begin() { return &(recv_buffer[0]); }
	char * recv_buffer_end() { return &(recv_buffer[0]) + bytes_received(); }

	// Buffer packing and unpacking
	template<class T> void pack(const int & send_to, const T &);
	template<class T> void unpack(char **, T &);

protected:
	// Hide copy constructor and assignment operator
	CommunicationBuffer(const CommunicationBuffer &);
	CommunicationBuffer & operator=(const CommunicationBuffer &);
	MPI_Request * send_requests;
	MPI_Request * recv_requests;
	std::vector<plint> proc_list;

private:
	virtual void before_send() = 0;
	virtual void after_receive() = 0;
	virtual void get_recv_size(std::vector<plint> & send_size, std::vector<plint>& recv_size) = 0;

	void allocate_request_arrays();
	void free_request_arrays();

	// MPI requests
	bool transfer_in_progress;
	int send_request_count;
	int recv_request_count;

	BufferType send_buffer;
	std::vector<char> recv_buffer;
};

class CommunicationBufferWithEqualSendAndReceiveSize : public CommunicationBuffer {
public:
	virtual void get_recv_size(std::vector<plint> &, std::vector<plint>&);
};

class CommunicationBufferCommunicatedSendAndReceiveSize : public CommunicationBuffer {
public:
	virtual void get_recv_size(std::vector<plint> &, std::vector<plint>&);
};

/*
 * Communicator for the Fsi force evaluation
 */
// Forward declaration
template<class T> class SolidNodeList;

template<class T>
struct VelocityInterpolationNode {
	Array<T, 3> delta_u;
	Array<T, 3> pos, pos_rel;
	plint obj_id;
};

template<class T, class Dirac, class Arithmetic>
class FsiForceCommunicator : public CommunicationBufferWithEqualSendAndReceiveSize {
public:
	typedef typename std::vector<VelocityInterpolationNode<T> > DataContainerType;
	typedef typename DataContainerType::iterator IteratorType;
	typedef std::vector<std::pair<geo::Rect<T>, plint> > DomainProcMapType;

	FsiForceCommunicator(const Arithmetic & a) : arithmetic(a) { }
	virtual ~FsiForceCommunicator() { }

	void set_domain_proc_map(const DomainProcMapType & map) { domain_proc_map = map; }

	void add_data(const VelocityInterpolationNode<T> & d) { data.push_back(d); }
	void clear_data() { data.clear(); }
	IteratorType data_begin() { return data.begin(); }
	IteratorType data_end() { return data.end(); }

private:
	virtual void before_send();
	virtual void after_receive();

	DataContainerType data;
	Arithmetic arithmetic;

	DomainProcMapType domain_proc_map;	// Neighboring processors and their corresponding domains
};

/*
 * Communicator for the particle force interpolation
 */
template<class T, class Arithmetic>
class ParticleForceCommunicator : public CommunicationBufferWithEqualSendAndReceiveSize {
public:
	enum ForceType { TmpForce, CollForce, FsiForce };

	typedef typename std::map<plint, RigidParticle3D<T> * > ObjMapType;
	typedef std::vector<std::pair<geo::Rect<T>, plint> > DomainProcMapType;

	virtual ~ParticleForceCommunicator() { }
	ParticleForceCommunicator(ObjMapType & obj_map, const Arithmetic & a) : obj_map(obj_map), arithmetic(a) { }

	void set_domain_proc_map(const DomainProcMapType & map) { domain_proc_map = map; }
	void pack_from_particles(ForceType);
	void unpack_to_particles(ForceType);

private:
	virtual void before_send() { }
	virtual void after_receive() { }

	Arithmetic arithmetic;
	DomainProcMapType domain_proc_map;	// Neighboring processors and their corresponding domains
	ObjMapType & obj_map;
};

/*
 * Communicator for particle state synchronization
 */
template<class T>
struct ParticleStateNode {
	Array<T, 3> pos, vel, ang_vel, force_coll, force_coll_last, torque_coll, torque_coll_last, force, torque, ang_momentum;
	Quaternion<T> orientation;
	T density, scale;
	plint shape_id, obj_id;
};

template<class T, class Arithmetic>
class ParticleStateCommunicator : public CommunicationBufferCommunicatedSendAndReceiveSize {
public:
	typedef std::vector<std::pair<geo::Rect<T>, plint> > DomainProcMapType;
	typedef typename std::vector<ParticleStateNode<T> > DataContainerType;
	typedef typename DataContainerType::iterator DataIterator;

	ParticleStateCommunicator(const Arithmetic & a) : arithmetic(a) { }

	void pack_particle(const RigidParticle3D<T> & particle);
	void set_domain_proc_map(const DomainProcMapType & d_p_map) { domain_proc_map = d_p_map; }

	void clear_data() { data.clear(); }
	DataIterator data_begin() { return data.begin(); }
	DataIterator data_end() { return data.end(); }

private:
	virtual void after_receive();
	virtual void before_send() { }

	Arithmetic arithmetic;
	DataContainerType data;
	DomainProcMapType domain_proc_map;
};

} /* namespace fsi */

} /* namespace plb */



#endif /* COMMUNICATIONBUFFER_H_ */
