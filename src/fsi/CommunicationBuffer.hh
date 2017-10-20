#ifndef COMMUNICATIONBUFFER_HH_
#define COMMUNICATIONBUFFER_HH_
#include "CommunicationBuffer.h"
#include "utils.h"
#include <stdexcept>

namespace plb {

namespace fsi {

/*** LazyValue ***/
template<class T>
void CommunicationBuffer::LazyValue<T>::finalize()
{
	PLB_PRECONDITION(vec != 0)
	utils::pack(&((*vec)[pos]), val);
}

/*** CommunicationBuffer ***/
void CommunicationBuffer::send_and_receive_no_wait(bool communicate_send_size)
{
	if(proc_list.empty())
		return;

	if(transfer_in_progress) {
		pcerr << "Attempting to send from an already sending CommunicationManager" << std::endl;
		exit(-1);
	}

	// Compute the size of the messages to send
	std::vector<plint> send_size(proc_list.size());

	for(plint i = 0; i < proc_list.size(); ++i) {
		BufferIterator it = send_buffer.find(proc_list[i]);
		if(it == send_buffer.end())
			send_size[i] = 0;
		else
			send_size[i] = it->second.size();
	}

	if(communicate_send_size) {
		// Send and receive sizes
		for(plint i = 0; i < proc_list.size(); ++i) {
			global::mpi().iSend(&(send_size[i]), 1, proc_list[i], &send_requests[i]);
			global::mpi().iRecv(&(recv_size_[i]), 1, proc_list[i], &recv_requests[i]);
		}

		// Wait for data transfer
		MPI_Waitall(recv_size_.size(), recv_requests, MPI_STATUSES_IGNORE);
		MPI_Waitall(send_size.size(), send_requests, MPI_STATUSES_IGNORE);
	}

	// Compute total data size to receive
	plint data_size = std::accumulate(recv_size_.begin(), recv_size_.end(), 0);

	// Reset the receive and send request counts
	recv_request_count = 0;
	send_request_count = 0;

	// Receive
	recv_buffer.resize(data_size);
	char * buffer_ptr = &(recv_buffer[0]);
	for(plint i = 0; i < proc_list.size(); ++i) {
		// Receive
		if(recv_size_[i] != 0) {
			//std::cout << global::mpi().getRank() << " <- " << proc_list[i] << " (" << recv_size_[i] << " bytes)"  << std::endl;
			global::mpi().iRecv(buffer_ptr,
								recv_size_[i],
								proc_list[i],
								&recv_requests[recv_request_count++]);
			buffer_ptr += recv_size_[i];
		}
	}

	// Send
	for(BufferIterator it = send_buffer.begin(); it != send_buffer.end(); ++it) {
		if(it->second.size() != 0) {
			//std::cout << global::mpi().getRank() << " -> " << it->first << " (" << it->second.size() << " bytes) " << std::endl;
			global::mpi().iSend(&(it->second[0]),
								it->second.size(),
								it->first,
								&send_requests[send_request_count++]);
		}
	}

	transfer_in_progress = true;
}

void CommunicationBuffer::finalize_send_and_receive()
{
	if(proc_list.empty())
		return;

	if( ! transfer_in_progress) {
		pcerr << "Attempting to finalize a send and receive from an idle CommunicationManager" << std::endl;
		exit(-1);
	}

	// Wait for all data transfer to finish
	//if(recv_request_count != 0)
	MPI_Waitall( recv_request_count, recv_requests, MPI_STATUSES_IGNORE );
	//if(send_request_count != 0)
	MPI_Waitall( send_request_count, send_requests, MPI_STATUSES_IGNORE );
	transfer_in_progress = false;
}

void CommunicationBuffer::clear_send_buffer()
{
	for(BufferIterator it = send_buffer.begin(); it != send_buffer.end(); ++it)
		it->second.clear();
}

void CommunicationBuffer::set_proc_list(const std::vector<plint> & proc_list)
{
	this->proc_list = proc_list;
	allocate_request_arrays();
}

void CommunicationBuffer::allocate_request_arrays()
{
	free_request_arrays();
	send_requests = new MPI_Request[proc_list.size()];
	recv_requests = new MPI_Request[proc_list.size()];
	recv_size_.resize(proc_list.size());
}

void CommunicationBuffer::free_request_arrays()
{
	if(send_requests)
		delete [] send_requests;
	if(recv_requests)
		delete [] recv_requests;
}

template<class T>
void CommunicationBuffer::pack(const int & send_to, const T & data)
{
	utils::pack(send_buffer[send_to], data);
}

template<class T>
void CommunicationBuffer::pack(const int & send_to, LazyValue<T> & data)
{
	plint insert_pos = send_buffer[send_to].size();

	// Pack a placeholder value (just so we have available memory later)
	this->pack(send_to, T());

	// Save the pointer to where the data began
	data.vec = &(send_buffer.at(send_to));
	data.pos = insert_pos;
}

template<class T>
void CommunicationBuffer::delegate_pack(const int & send_to, const T & packer)
{
	packer.pack(send_buffer.at(send_to));
}

template<class T>
void CommunicationBuffer::unpack(char *& buffer, T & data)
{
	utils::unpack(buffer, data);
}

void CommunicationBuffer::set_recv_size(plint i, plint recv_size)
{
	std::vector<plint>::iterator it = std::find(proc_list.begin(), proc_list.end(), i);
	PLB_PRECONDITION(it != proc_list.end())
	plint ind = it - proc_list.begin();
	recv_size_[ind] = recv_size;
}

} /* namespace fsi */

} /* namespace plb */


#endif /* COMMUNICATIONBUFFER_HH_ */
