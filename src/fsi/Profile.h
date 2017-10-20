/*
 * Profile.h
 *
 *  Created on: Apr 16, 2014
 *      Author: niber
 */

#ifndef PROFILE_H_
#define PROFILE_H_
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "core/globalDefs.h"

namespace plb {

namespace fsi {

struct ProfileEntry {
	double start_time;
	double total_time;
	pluint invoke_count;
	bool is_running;
};

class Profile {
public:
	typedef std::map<std::string, ProfileEntry> ContainerType;
	typedef ContainerType::iterator IteratorType;

	static void start_timer(const std::string & name)
	{
#ifdef FSI_PROFILE
		IteratorType it = entries.find(name);
		// Create new entry if none exists with the given name
		if(it == entries.end()) {
			ProfileEntry entry;
			entry.start_time = 0;
			entry.total_time = 0;
			entry.invoke_count = 0;
			entry.is_running = false;

			std::pair<IteratorType, bool> res =
					entries.insert(std::make_pair(name, entry));
			it = res.first;
		}

		PLB_PRECONDITION( ! it->second.is_running)
		it->second.is_running = true;
		it->second.invoke_count += 1;
		it->second.start_time = MPI_Wtime();
#endif
	}

	static void stop_timer(const std::string & name)
	{
#ifdef FSI_PROFILE
		std::map<std::string, ProfileEntry>::iterator it = entries.find(name);
		if(it == entries.end())
			std::cerr << "Trying to stop a timer, " << name << ", that did not exist" << std::endl;

		PLB_PRECONDITION(it->second.is_running)
		it->second.is_running = false;
		it->second.total_time += MPI_Wtime() - it->second.start_time;
#endif
	}

	static void write_report(const std::string & file_name)
	{
#ifdef FSI_PROFILE
		std::stringstream fname;
		fname << file_name;
		fname << global::mpi().getRank();

		FileName name;
		name.setName(fname.str());
		name.setPath(global::directories().getOutputDir());
		name.setExt("txt");

		std::fstream out(name.get().c_str(), std::ios::out);
		for(IteratorType it = entries.begin(); it != entries.end(); ++it) {
			out << it->first << ": " << std::endl
				<< "  Number of calls: " << it->second.invoke_count << std::endl
				<< "  Total time: " << it->second.total_time << std::endl
				<< "  Average time: " << it->second.total_time / (double) it->second.invoke_count << std::endl;
		}
		out.close();
#endif
	}

private:
	Profile() { }
	~Profile() { }
#ifdef FSI_PROFILE
	static std::map<std::string, ProfileEntry> entries;
#endif
};

#ifdef FSI_PROFILE
std::map<std::string, ProfileEntry> Profile::entries;
#endif

}

}

#endif /* PROFILE_H_ */
