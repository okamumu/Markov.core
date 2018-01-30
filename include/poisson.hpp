/*
  poisson
 */

namespace marlib {

	template <typename VectorT>
	class poisson {
		static constexpr double NORMALQ_LOWER_Q = 3.0;
		static constexpr double NORMALQ_UPPER_Q = 37.0;

		static constexpr double NORMALQ_LOWER_LOGP = -689.0;
		static constexpr double NORMALQ_UPPER_LOGP = -6.6;

		static constexpr double NORMALQ_EPS = 1.0e-8;

		static constexpr double LOG2PIOVER2 = 0.9189385332046727417803297364; // log(2pi) / 2

		static constexpr double POISSON_LAMBDA_MIN = 3.0;
		static constexpr int POISSON_RIGHT_MAX = 23;

	public:
    poisson(int size, const VectorT& prob);
		~poisson();

		int left() const;
		int right() const;
		double weight() const;
    double lambda() const;

		double operator()(int i) const;

    static int rightbound(double lambda, double eps);

    void set(double lambda, int left, int right);

	private:
		double m_lambda;
		int m_size;
		int m_left;
		int m_right;
		VectorT m_prob;
		double m_weight;

		static double normalt(double x);
		static double normalq(double p);
		double set_prob();

	};

	template <typename VectorT>
	poisson<VectorT>::poisson(int size, const VectorT& prob)
	: m_size(size), m_left(0), m_right(size-1), m_prob(prob) {}

	template <typename VectorT>
	poisson<VectorT>::~poisson() {}

  template <typename VectorT>
	inline int poisson<VectorT>::left() const {
		return m_left;
	}

	template <typename VectorT>
	inline int poisson<VectorT>::right() const {
		return m_right;
	}

	template <typename VectorT>
	inline double poisson<VectorT>::weight() const {
		return m_weight;
	}

	template <typename VectorT>
	inline double poisson<VectorT>::lambda() const {
		return m_lambda;
	}

	template <typename VectorT>
	inline double poisson<VectorT>::operator()(int i) const {
		return m_prob[i - m_left];
	}

	template <typename VectorT>
	void poisson<VectorT>::set(double lambda, int left, int right) {
		assert(m_size >= right - left + 1);
		m_lambda = lambda;
		m_left = left;
		m_right = right;
		m_weight = set_prob();
	}

	/*
	Description
	return a tail probability of standard normal distribution
	Parameters
	IN
	x: input value
	OUT
	return value
	*/

	template <typename VectorT>
	double poisson<VectorT>::normalt(double x) {
		double x2 = x*x;
		double tmp = x;
		double sum = 1 / tmp;
		tmp = tmp * x2;
		sum = sum - 1 / tmp;
		tmp = tmp * x2;
		sum = sum + 3 / tmp;
		tmp = tmp * x2;
		sum = sum - 15 / tmp;
		tmp = tmp * x2;
		sum = sum + 105 / tmp;
		return (log(sum) - x2/2 - LOG2PIOVER2);
	}

	/*
	Description
	return a quantile of standard normal distribution
	Parameters
	IN
	p: probability for the quantile
	OUT
	return value
	*/

	template <typename VectorT>
	double poisson<VectorT>::normalq(double p) {
		const double leps = log(p);
		assert(leps <= NORMALQ_UPPER_LOGP && leps >= NORMALQ_LOWER_LOGP);
		double l = NORMALQ_LOWER_Q;
		double u = NORMALQ_UPPER_Q;
		double m = (l + u) / 2;
		double fm = normalt(m) - leps;
		while (std::fabs(fm) > NORMALQ_EPS) {
			if (fm > 0) {
				l = m;
			} else {
				u = m;
			}
			m = (l + u)/2;
			fm = normalt(m) - leps;
		}
		return m;
	}

	/*
	! Description: compute the right bound of Poisson range
	!              for a given error tolerance
	!
	! Parameters:
	!   IN
	!    lambda: Poisson rate (mean)
	!    epsi: error tolerance
	!   OUT
	!    right bound is a return value
	*/

	template <typename VectorT>
	int poisson<VectorT>::rightbound(double lambda, double eps) {
		int right;
		if (lambda == 0) {
			return 0;
		}
		if (lambda < POISSON_LAMBDA_MIN) {
			double tmp = exp(-lambda);
			double total = tmp;
			right = 0;
			for (int k=1; k<=POISSON_RIGHT_MAX; k++) {
				right++;
				tmp *= lambda / right;
				total += tmp;
				if (total + eps >= 1)
				break;
			}
		} else {
			double z = normalq(eps);
			double tmp = z + sqrt(4 * lambda - 1);
			right = static_cast<int>(tmp * tmp / 4 + 1);
		}
		return right;
	}

	/*
	! Description & parameters :
	!  IN
	!    lambda: Poisson parameter (mean)
	!    left, right: left and right bounds
	!  OUT
	!    prob: Poisson probabilities from left to right bounds
	!    weight: weight value, i.e., exact pmf is prob[i]/weight
	*/

	template <typename VectorT>
	double poisson<VectorT>::set_prob() {
		int mode = static_cast<int>(m_lambda);
		if (mode >= 1) {
			m_prob[mode-m_left] = exp(-m_lambda + mode * log(m_lambda) - LOG2PIOVER2 - (mode + 1.0/2.0) * log(mode) + mode);
		} else {
			m_prob[mode-m_left] = exp(-m_lambda);
		}
		// -- down --
		for (int j=mode; j>=m_left+1; j--) {
			m_prob[j-1-m_left] = m_prob[j-m_left] * j / m_lambda;
		}
		// -- up --
		for (int j=mode; j<=m_right-1; j++) {
			m_prob[j+1-m_left] = m_prob[j-m_left] * m_lambda / (j+1);
		}
		// -- compute W --
		double weight = 0;
		int s = m_left;
		int t = m_right;
		while (s < t) {
			if (m_prob[s-m_left] <= m_prob[t-m_left]) {
				weight += m_prob[s-m_left];
				s++;
			} else {
				weight += m_prob[t-m_left];
				t--;
			}
		}
		weight += m_prob[s-m_left];
		return weight;
	}
}
