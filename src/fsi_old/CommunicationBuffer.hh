#ifndef COMMUNICATIONBUFFER_HH_
#define COMMUNICATIONBUFFER_HH_
#include "CommunicationBuffer.h"
#include <stdexcept>

namespace plb {

namespace fsi {

void CommunicationBuffer::send_and_receive_no_wait()
{
	// Call the before send callback (implemented in the subclasses)
	before_send();

	if(proc_list.empty())
		return;

	if(transfer_in_progress) {
		pcerr << "Attempting to send from an already sending CommunicationManager" << std::endl;
		exit(-1);
	}

	// Compute the size of the messages to send
	std::vector<plint> send_size(proc_list.size());
	std::vector<plint> recv_size(proc_list.size());

	for(plint i = 0; i < proc_list.size(); ++i) {
		BufferIterator it = send_buffer.find(proc_list[i]);
		if(it == send_buffer.end())
			send_size[i] = 0;
		else
			send_size[i] = it->second.size();
	}

	// Get the size to receive from each of the other processors
	get_recv_size(send_size, recv_size);

	// Compute total data size to receive
	plint data_size = std::accumulate(recv_size.begin(), recv_size.end(), 0);

	// Reset the receive and send request counts
	recv_request_count = 0;
	send_request_count = 0;

	// Receive
	recv_buffer.resize(data_size);
	char * buffer_ptr = &(recv_buffer[0]);
	for(plint i = 0; i < proc_list.size(); ++i) {
		// Receive
		if(recv_size[i] != 0) {
			/*std::cout << global::mpi().getRank() << " <- " << proc_list[i]
						<< " (" << recv_size[i] << " bytes)"  << std::endl;*/
			global::mpi().iRecv(buffer_ptr,
								recv_size[i],
								proc_list[i],
								&recv_requests[recv_request_count++]);
			buffer_ptr += recv_size[i];
		}
	}

	// Send
	for(BufferIterator it = send_buffer.begin(); it != send_buffer.end(); ++it) {
		if(it->second.size() != 0) {
			/*std::cout << global::mpi().getRank() << " -> " << it->first
						<< " (" << it->second.size() << " bytes) " << std::endl;*/
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

	after_receive();
}

void CommunicationBuffer::clear_send_buffer()
{
	for(BufferIterator it = send_buffer.begin(); it != send_buffer.end(); ++it)
		it->second.clear();
}

void CommunicationBuffer::allocate_request_arrays()
{
	free_request_arrays();
	send_requests = new MPI_Request[proc_list.size()];
	recv_requests = new MPI_Request[proc_list.size()];
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
	const char * const buf = reinterpret_cast<const char *>(&data);
	send_buffer[send_to].insert(send_buffer[send_to].end(), buf, buf + sizeof(T));
}

template<class T>
void CommunicationBuffer::unpack(char ** buffer, T & data)
{
	std::copy(*buffer, (*buffer)+sizeof(T), reinterpret_cast<char *>(&data));
	*buffer += sizeof(T);
}

/*
 * CommunicationBufferWithEqualSendAndReceiveSize
 */
void CommunicationBufferWithEqualSendAndReceiveSize::get_recv_size(
		std::vector<plint> & send_size,
		std::vector<plint> & recv_size)
{
	for(pluint i = 0; i < send_size.size(); ++i)
		recv_size[i] = send_size[i];
}

/*
 * CommunicationBufferCommunicatedSendAndReceiveSize
 */
void CommunicationBufferCommunicatedSendAndReceiveSize::get_recv_size(
		std::vector<plint> & send_size,
		std::vector<plint> & recv_size)
{
	// Send and receive sizes
	for(plint i = 0; i < proc_list.size(); ++i) {
		global::mpi().iSend(&(send_size[i]), 1, proc_list[i], &send_requests[i]);
		global::mpi().iRecv(&(recv_size[i]), 1, proc_list[i], &recv_requests[i]);
	}

	// Wait for data transfer
	MPI_Waitall(recv_size.size(), recv_requests, MPI_STATUSES_IGNORE);
	MPI_Waitall(send_size.size(), send_requests, MPI_STATUSES_IGNORE);
}

/*
 * FsiForceCommunicator
 */
template<class T, class Dirac, class Arithmetic>
void FsiForceCommunicator<T, Dirac, Arithmetic>::before_send()
{
	Box3D d_box, p_box;

	clear_send_buffer();

	// Pack all data into the send buffer
	for(pluint i = 0, i_end = data.size(); i < i_end; ++i) {
		const VelocityInterpolationNode<T> & node = data[i];

		// Determine which processors that should receive the data
		typename std::vector<std::pair<geo::Rect<T>, plint> >::iterator it = domain_proc_map.begin(),
				it_end = domain_proc_map.end();
		for(; it != it_end; ++it) {
			// Is the compact support of this node in the currently probed processors domain?
			if(it->first.contains_or_on_boundary(node.pos, arithmetic)) {
				// Then pack the information into the send buffer for that processor
				pack(it->second, node.obj_id);
				pack(it->second, node.pos);
				pack(it->second, node.pos_rel);
				pack(it->second, node.delta_u);
			}
		}
	}
}

template<class T, class Dirac, class Arithmetic>
void FsiForceCommunicator<T, Dirac, Arithmetic>::after_receive()
{
	VelocityInterpolationNode<T> node;
	clear_data();

	for(char * buffer_ptr = recv_buffer_begin(), * buffer_end = recv_buffer_end(); buffer_ptr != buffer_end; ) {
		unpack(&buffer_ptr, node.obj_id);
		unpack(&buffer_ptr, node.pos);
		unpack(&buffer_ptr, node.pos_rel);
		unpack(&buffer_ptr, node.delta_u);
		data.push_back(node);
	}
}

/*
 * ParticleForceCommunicator
 */
template<class T, class Arithmetic>
void ParticleForceCommunicator<T, Arithmetic>::pack_from_particles(ForceType force_type)
{
	clear_send_buffer();

	RigidParticle3D<T> * particle;
	geo::Sphere<T> b_sphere;

	// The switch statement is put outside of the inner loop to avoid unnecessary branching
	switch(force_type) {
	case TmpForce:
		// Pack force and torque data
		for(typename ObjMapType::iterator it = obj_map.begin(); it != obj_map.end(); ++it) {
			particle = it->second;
			b_sphere = particle->get_bounding_sphere();

			for(typename DomainProcMapType::iterator dom_it = domain_proc_map.begin(); dom_it != domain_proc_map.end(); ++dom_it) {
				if(geo::does_intersect(b_sphere, dom_it->first, arithmetic)) {
					pack(dom_it->second, it->first);
					pack(dom_it->second, particle->fsi_force_tmp);
					pack(dom_it->second, particle->fsi_torque_tmp);
				}
			}
		}
		break;
	case CollForce:
		// Pack force and torque data
		for(typename ObjMapType::iterator it = obj_map.begin(); it != obj_map.end(); ++it) {
			particle = it->second;
			b_sphere = particle->get_bounding_sphere();

			for(typename DomainProcMapType::iterator dom_it = domain_proc_map.begin(); dom_it != domain_proc_map.end(); ++dom_it) {
				if(geo::does_intersect(b_sphere, dom_it->first, arithmetic)) {
					pack(dom_it->second, it->first);
					pack(dom_it->second, particle->coll_force);
					pack(dom_it->second, particle->coll_torque);
				}
			}
		}
		break;
	}
}

template<class T, class Arithmetic>
void ParticleForceCommunicator<T, Arithmetic>::unpack_to_particles(ForceType force_type)
{
	Array<T, 3> force, torque;
	typename ObjMapType::key_type id;
	RigidParticle3D<T> * particle;

	// Unpack and add the contributions
	switch(force_type) {
	case TmpForce:
		for(char * buffer_ptr = recv_buffer_begin(), * buffer_end = recv_buffer_end(); buffer_ptr != buffer_end; ) {
			unpack(&buffer_ptr, id);
			unpack(&buffer_ptr, force);
			unpack(&buffer_ptr, torque);

			PLB_PRECONDITION(obj_map.find(id) != obj_map.end());

			particle = obj_map.at(id);
			particle->fsi_force_tmp += force;
			particle->fsi_torque_tmp += torque;
		}
		break;
	case CollForce:
		for(char * buffer_ptr = recv_buffer_begin(), * buffer_end = recv_buffer_end(); buffer_ptr != buffer_end; ) {
			unpack(&buffer_ptr, id);
			unpack(&buffer_ptr, force);
			unpack(&buffer_ptr, torque);

			PLB_PRECONDITION(obj_map.find(id) != obj_map.end());

			particle = obj_map.at(id);
			particle->coll_force += force;
			particle->coll_torque += torque;
		}
		break;
	}
}

/*
 * ParticleStateCommunicator
 */
template<class T, class Arithmetic>
void ParticleStateCommunicator<T, Arithmetic>::pack_particle(const RigidParticle3D<T> & particle)
{
	// Particle bounding sphere
	geo::Sphere<T> b_sphere = particle.get_bounding_sphere();

	// Pack to the appropriate processors
	for(typename DomainProcMapType::iterator it = domain_proc_map.begin(); it != domain_proc_map.end(); ++it) {
		if(geo::does_intersect(b_sphere, it->first, arithmetic)) {
			pack(it->second, particle.get_id());
			pack(it->second, particle.get_shape_id());
			pack(it->second, particle.get_angular_velocity());
			pack(it->second, particle.get_orientation());
			pack(it->second, particle.get_position());
			pack(it->second, particle.get_velocity());
			pack(it->second, particle.get_density());
			pack(it->second, particle.get_scale());
			pack(it->second, particle.coll_force);
			pack(it->second, particle.coll_force_last);
			pack(it->second, particle.coll_torque);
			pack(it->second, particle.coll_torque_last);
			pack(it->second, particle.get_force());
			pack(it->second, particle.get_torque());
			pack(it->second, particle.ang_momentum);
		}
	}
}

template<class T, class Arithmetic>
void ParticleStateCommunicator<T, Arithmetic>::after_receive()
{
	// Unpack data
	data.clear();
	ParticleStateNode<T> node;

	for(char * buffer_ptr = recv_buffer_begin(), * buffer_end = recv_buffer_end(); buffer_ptr != buffer_end; ) {
		unpack(&buffer_ptr, node.obj_id);
		unpack(&buffer_ptr, node.shape_id);
		unpack(&buffer_ptr, node.ang_vel);
		unpack(&buffer_ptr, node.orientation);
		unpack(&buffer_ptr, node.pos);
		unpack(&buffer_ptr, node.vel);
		unpack(&buffer_ptr, node.density);
		unpack(&buffer_ptr, node.scale);
		unpack(&buffer_ptr, node.force_coll);
		unpack(&buffer_ptr, node.force_coll_last);
		unpack(&buffer_ptr, node.torque_coll);
		unpack(&buffer_ptr, node.torque_coll_last);
		unpack(&buffer_ptr, node.force);
		unpack(&buffer_ptr, node.torque);
		unpack(&buffer_ptr, node.ang_momentum);
		data.push_back(node);
	}
}

} /* namespace fsi */

} /* namespace plb */


#endif /* COMMUNICATIONBUFFER_HH_ */
