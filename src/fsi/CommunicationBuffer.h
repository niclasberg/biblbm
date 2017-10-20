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

	template<class T>
	class LazyValue {
	public:
		LazyValue() : val(), vec(0), pos(0) { }
		LazyValue(T v) : val(v), vec(0), pos(0) { }

		T & get() { return val; }
		const T & get() const { return val; }

		void finalize();

	private:
		T val;
		std::vector<char> * vec;
		plint pos;
		friend class CommunicationBuffer;
	};

	CommunicationBuffer()
	: transfer_in_progress(false), proc_list(), send_requests(0), recv_requests(0),
	  recv_request_count(0), send_request_count(0), recv_size_()
	{

	}

	CommunicationBuffer(const std::vector<plint> & proc_list)
	: transfer_in_progress(false), proc_list(proc_list), send_requests(0), recv_requests(0)
	{
		allocate_request_arrays();
	}

	virtual ~CommunicationBuffer() { free_request_arrays(); }

	void set_proc_list(const std::vector<plint> & proc_list);

	void send_and_receive_no_wait(bool communicate_send_size);
	void finalize_send_and_receive();

	void clear_send_buffer();

	pluint bytes_received() const { return recv_buffer.size(); }
	char * recv_buffer_begin() { return &(recv_buffer[0]); }
	char * recv_buffer_end() { return &(recv_buffer[0]) + bytes_received(); }

	// Buffer packing and unpacking
	template<class T> void pack(const int & send_to, const T &);
	template<class T> void pack(const int & send_to, LazyValue<T> &);
	template<class T> void delegate_pack(const int & send_to, const T &);
	template<class T> void unpack(char *&, T &);

	// Receive sizes
	void set_recv_size(plint i, plint recv_size);

protected:
	// Hide copy constructor and assignment operator
	CommunicationBuffer(const CommunicationBuffer &);
	CommunicationBuffer & operator=(const CommunicationBuffer &);
	MPI_Request * send_requests;
	MPI_Request * recv_requests;
	std::vector<plint> proc_list;
	std::vector<plint> recv_size_;

private:

	void allocate_request_arrays();
	void free_request_arrays();

	// MPI requests
	bool transfer_in_progress;
	int send_request_count;
	int recv_request_count;

	BufferType send_buffer;
	std::vector<char> recv_buffer;
};

} /* namespace fsi */

} /* namespace plb */



#endif /* COMMUNICATIONBUFFER_H_ */
