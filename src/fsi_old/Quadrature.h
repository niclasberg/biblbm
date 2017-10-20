#ifndef QUADRATURE_H_
#define QUADRATURE_H_

namespace plb {

namespace fsi {

template<class T>
class Quadrature {
public:
	Quadrature() { set_order(1); }
	Quadrature(plint order_) { set_order(order_); }

	void set_order(plint);
	const Array<T, 3> & get_point(plint i) const { return points[i]; }
	const T & get_weight(plint i) const { return weights[i]; }
	plint get_num_points() const { return points.size(); }
private:
	plint order;
	std::vector<T> weights;
	std::vector<Array<T, 2> > points;
};

template<class T>
inline void Quadrature<T>::set_order(plint order_)
{
	weights.clear();
	points.clear();

	switch(order_) {
	case 1:
		weights.push_back((T)0.5);
		points.push_back(Array<T, 2>((T)1.0/3.0, (T)1.0/3.0));
		break;
	case 2:
		weights.push_back(1.0/6.0);
		points.push_back(Array<T, 2>((T)1.0/6.0, (T)1.0/6.0));

		weights.push_back(1.0/6.0);
		points.push_back(Array<T, 2>((T)2.0/3.0, (T)1.0/6.0));

		weights.push_back(1.0/6.0);
		points.push_back(Array<T, 2>((T)1.0/6.0, (T)2.0/3.0));
		break;
	case 3:
		weights.push_back(-9.0/32.0);
		points.push_back(Array<T, 2>((T)1.0/3.0, (T)1.0/3.0));

		weights.push_back(25.0/96.0);
		points.push_back(Array<T, 2>((T)3.0/5.0, (T)1.0/5.0));

		weights.push_back(25.0/96.0);
		points.push_back(Array<T, 2>((T)1.0/5.0, (T)3.0/5.0));

		weights.push_back(25.0/96.0);
		points.push_back(Array<T, 2>((T)1.0/5.0, (T)1.0/5.0));
		break;
	case 4:
		weights.push_back(1.0/40.0);
		points.push_back(Array<T, 2>((T)0, (T)0));

		weights.push_back(1.0/15.0);
		points.push_back(Array<T, 2>((T)1.0/2.0, (T)0));

		weights.push_back(1.0/40.0);
		points.push_back(Array<T, 2>((T)1, (T)0));

		weights.push_back(1.0/15.0);
		points.push_back(Array<T, 2>((T)1./2., (T)1./2.));

		weights.push_back(1.0/40.0);
		points.push_back(Array<T, 2>((T)0, (T)1));

		weights.push_back(1.0/15.0);
		points.push_back(Array<T, 2>((T)0, (T)1./2.));

		weights.push_back(9.0/40.0);
		points.push_back(Array<T, 2>((T)1./3., (T)1./3.));
		break;
	default:
		pcerr << "Unknown quadrature order " << order_ << std::endl;
		return;
	}
	order = order_;
}

}

}




#endif /* QUADRATURE_H_ */
