/*
 * ZBuffer.h
 *
 *  Created on: Jul 23, 2015
 *      Author: niber
 */

#ifndef ZBUFFER_H_
#define ZBUFFER_H_
#include <vector>
#include "core/geometry3D.h"
#include "utils.h"

namespace plb {

namespace fsi {

template<class T, class Arithmetic>
class ZBuffer {
public:
	typedef typename std::vector<T>::const_iterator const_iterator;

	ZBuffer(const Arithmetic & arithmetic)
	: x0_(0), y0_(0), nx_(0), ny_(0), data_(), arithmetic_(arithmetic)
	 { }

	ZBuffer(const Arithmetic & arithmetic, plint x0, plint x1, plint y0, plint y1)
	: x0_(x0), y0_(y0), nx_(x1-x0+1), ny_(y1-y0+1), data_(), arithmetic_(arithmetic)
	{
		data_.resize(count_cells());
	}

	void clear()
	{
		for(plint i = 0; i < data_.size(); ++i)
			data_[i].clear();
	}

	void put(plint iX, plint iY, const T & val)
	{
		PLB_PRECONDITION(idx(iX, iY) >= 0 && idx(iX, iY) < count_cells())
		T z = val;
		arithmetic_.remap_position_z(z);
		if(idx(iX, iY) < 0 || idx(iX, iY) >= count_cells())
			std::cerr << "ERROR: ZBuffer::put: Attempting to insert at (x, y) =  (" << iX << ", " << iY << ")" <<
				" local coordinates: " << idx_x(iX) << ", " << idx_y(iY) <<
				" (domain (x0, y0) = (" << x0_ << ", " << y0_ << "), nx = " << nx_ << " ny = " << ny_ << std::endl;
		data_[idx(iX, iY)].push_back(z);
	}

	std::pair<const_iterator, const_iterator> get(plint iX, plint iY) const
	{
		plint ind = idx(iX, iY);
		return std::make_pair(data_[ind].begin(), data_[ind].end());
	}

	void sort_entries()
	{
		for(plint i = 0; i < data_.size(); ++i)
			std::sort(data_[i].begin(), data_[i].end());
	}

	template<class Compare>
	void sort_entries(const Compare & comp)
	{
		for(plint i = 0; i < data_.size(); ++i)
			std::sort(data_[i].begin(), data_[i].end(), comp);
	}

	void set_extent(plint x0, plint x1, plint y0, plint y1)
	{
		nx_ = std::min(x1 - x0 + 1, arithmetic_.get_nx());
		ny_ = std::min(y1 - y0 + 1, arithmetic_.get_ny());

		set_location(x0, y0);
		PLB_PRECONDITION(count_cells() >= 0)

		data_.resize(count_cells());
	}

	void set_location(plint x0, plint y0)
	{
		x0_ = x0;
		y0_ = y0;

		//arithmetic_.remap_index_x(x0_);
		//arithmetic_.remap_index_y(y0_);
	}

	plint count() const
	{
		plint ret = 0;
		for(plint i = 0; i < nx()*ny(); ++i)
			ret += data_[i].size();
		return ret;
	}

	/*
	 * Demands: x1 > x0 and y1 > y0
	 */
	template<class BufferType>
	void pack(plint x0, plint x1, plint y0, plint y1, BufferType & buf) const
	{
		// Get relative indices
		plint iX0 = idx_x(x0);
		plint iX1 = x1 + iX0 - x0;
		plint iY0 = idx_y(y0);
		plint iY1 = y1 + iY0 - y0;

		// Count number of entries to send
		plint num_entries = 0;
		plint iPmin = arithmetic_.periodic_x() ? -1 : 0;
		plint iPmax = arithmetic_.periodic_x() ? 1 : 0;
		plint jPmin = arithmetic_.periodic_y() ? -1 : 0;
		plint jPmax = arithmetic_.periodic_y() ? 1 : 0;
		plint nx_g = arithmetic_.get_nx();
		plint ny_g = arithmetic_.get_ny();

		for(plint i = iPmin; i <= iPmax; ++i)
			for(plint iX = std::max(iX0+i*nx_g, (plint)0); iX <= std::min(iX1 + i*nx_g, nx()-1); ++iX)
				for(plint j = jPmin; j <= jPmax; ++j)
					for(plint iY = std::max(iY0+j*ny_g, (plint)0); iY <= std::min(iY1 + j*ny_g, ny()-1); ++iY)
						num_entries += data_[iX*ny()+ iY].size();
		utils::pack(buf, num_entries);

		for(plint i = iPmin; i <= iPmax; ++i)
			for(plint iX = std::max(iX0+i*nx_g, (plint)0); iX <= std::min(iX1 + i*nx_g, nx()-1); ++iX) {
				plint iX2 = iX + x0_;
				arithmetic_.remap_index_x(iX2);

				for(plint j = jPmin; j <= jPmax; ++j)
					for(plint iY = std::max(iY0+j*ny_g, (plint)0); iY <= std::min(iY1 + j*ny_g, ny()-1); ++iY) {
						plint iY2 = iY + y0_;
						arithmetic_.remap_index_y(iY2);

						const std::vector<T> & d = data_[iX*ny()+ iY];

						plint size = d.size();

						if(size > 0) {
							utils::pack(buf, iX2);
							utils::pack(buf, iY2);
							utils::pack(buf, size);

							for(const_iterator it = d.begin(); it != d.end(); ++it)
								utils::pack(buf, *it);
						}
					}
			}
	}

	template<class BufferType>
	void unpack(BufferType & buf)
	{
		plint num_entries, iX, iY, size;
		T entry;
		plint unpacked = 0;
		utils::unpack(buf, num_entries);

		while(unpacked < num_entries) {
			utils::unpack(buf, iX);
			utils::unpack(buf, iY);
			utils::unpack(buf, size);

			//std::cout << "[" << x0_ << ", " << x1_ << "] " << iX << std::endl;

			std::vector<T> & d = data_[idx(iX, iY)];
			for(plint i = 0; i < size; ++i) {
				utils::unpack(buf, entry);
				d.push_back(entry);
			}
			unpacked += size;
		}
	}

	void dump(std::ostream & out) const
	{
		for(plint iX = 0; iX < nx_; ++iX) {
			plint iX2 = iX + x0_;
			arithmetic_.remap_index_x(iX2);

			for(plint iY = 0; iY < ny_; ++iY) {
				plint iY2 = iY + y0_;
				arithmetic_.remap_index_y(iY2);

				const std::vector<T> & d = data_[iX*ny()+ iY];
				plint size = d.size();

				if(size > 0) {
					out << iX2 << " " << iY2 << " " << size << " ";

					for(const_iterator it = d.begin(); it != d.end(); ++it)
						out << *it << " ";
				}
			}
		}
	}

	plint nx() const { return nx_; }
	plint ny() const { return ny_; }
	plint count_cells() const { return nx()*ny(); }

	plint x0() const { return x0_; }
	plint y0() const { return y0_; }
	plint x1() const { return x0_ + nx_ - 1; }
	plint y1() const { return y0_ + ny_ - 1; }

private:
	plint idx_x(plint iX) const { return arithmetic_.int_diff_x(iX, x0_); }
	plint idx_y(plint iY) const { return arithmetic_.int_diff_y(iY, y0_); }
	plint idx(plint iX, plint iY) const { return idx_x(iX)*ny() + idx_y(iY); }

private:
	plint x0_, y0_;
	plint nx_, ny_;
	Arithmetic arithmetic_;
	std::vector<std::vector<T> > data_;
};

}

}



#endif /* ZBUFFER_H_ */
