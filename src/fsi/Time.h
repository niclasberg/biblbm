/*
 * Time.h
 *
 *  Created on: Apr 20, 2015
 *      Author: niber
 */

#ifndef TIME_H_
#define TIME_H_

namespace plb {

namespace fsi {

class Time {
public:
	Time() : it_(0) { }
	Time(plint it) : it_(it) { }

	void write_checkpoint(FileName file_name) const
	{
		file_name.defaultPath(global::directories().getOutputDir());
		std::ofstream out(file_name.get().c_str(), std::ios::out);
		out << it_;
		out.close();
	}

	void read_checkpoint(FileName file_name)
	{
		file_name.defaultPath(global::directories().getInputDir());
		std::ifstream in(file_name.get().c_str(), std::ios::in);
		in >> it_;
		in.close();
	}

	bool is_multiple_of(plint N) const { return (it_ % N) == 0; }
	bool operator<(plint it) const { return it_ < it; }
	bool operator<=(plint it) const { return it_ <= it; }
	bool operator==(plint it) const { return it_ == it; }
	bool operator!=(plint it) const { return it_ != it; }
	Time & operator++() { ++it_; return *this; }
	Time operator++(int) { Time old = *this; ++it_; return old; }

	plint to_plint() const { return it_; }

private:
	plint it_;
};

std::ostream & operator<<(std::ostream & ss, const Time & t)
{
	ss << t.to_plint();
	return ss;
}

}

}


#endif /* TIME_H_ */
