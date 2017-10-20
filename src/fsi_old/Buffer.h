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

	template<class T2>
	void pack(const T2 & val)
	{
		*((T2 *) data_ptr) = val;
		data_ptr += sizeof(T2);
	}

	template<class T2>
	void unpack(T2 & val)
	{
		val = *((T2 *) data_ptr);
		data_ptr += sizeof(T2);
	}

private:
	T * data_ptr;
	std::vector<T> data;
};


}

}




#endif /* BUFFER_H_ */
