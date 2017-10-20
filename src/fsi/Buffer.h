#ifndef BUFFER_H_
#define BUFFER_H_

namespace plb {

namespace fsi {

template<class T>
class Buffer {
public:
	Buffer() { rewind_ptr(); }
	Buffer(pluint size) : data(size) { rewind_ptr(); }
	void rewind_ptr() { data_ptr = &(data[0]); }
	void advance_ptr(plint bytes) { data_ptr += bytes; }
	T * get_data_ptr() { return data_ptr; }
	bool has_data() const { return (data_ptr - &(data[0])) < data.size(); }

	template<class T2>
	void pack(const T2 & val)
	{
		*((T2 *) data_ptr) = val;
		data_ptr += sizeof(T2)/sizeof(T);
	}

	template<class T2>
	void unpack(T2 & val)
	{
		val = *((T2 *) data_ptr);
		data_ptr += sizeof(T2)/sizeof(T);
	}

	plint size() const { return data.size(); }

private:
	T * data_ptr;
	std::vector<T> data;
};

class DynamicBuffer {
public:
	DynamicBuffer() { rewind_ptr(); }
	void rewind_ptr() { data_ptr = &(data[0]); }
	void advance_ptr(plint bytes) { data_ptr += bytes; }
	char * get_data_ptr() { return data_ptr; }

	template<class T2>
	void pack(T2 val)
	{
		char * v = reinterpret_cast<char *>(&val);
		data.insert(data.end(), v, v + sizeof(T2));
		advance_ptr(sizeof(T2));
	}

	template<class T2>
	void unpack(T2 & val)
	{
		val = *((T2 *) data_ptr);
		advance_ptr(sizeof(T2));
	}

	plint size() const { return data.size(); }

private:
	char * data_ptr;
	std::vector<char> data;
};

}

}


#endif /* BUFFER_H_ */
