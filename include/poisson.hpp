/*
  poisson
 */

namespace marlib {

	template <typename ValueT>
	class poisson {
		static constexpr ValueT NORMALQ_LOWER_Q = 3.0;
		static constexpr ValueT NORMALQ_UPPER_Q = 37.0;

		static constexpr ValueT NORMALQ_LOWER_LOGP = -689.0;
		static constexpr ValueT NORMALQ_UPPER_LOGP = -6.6;

		static constexpr ValueT NORMALQ_EPS = 1.0e-8;

		static constexpr ValueT LOG2PIOVER2 = 0.9189385332046727417803297364; // log(2pi) / 2

		static constexpr ValueT POISSON_LAMBDA_MIN = 3.0;
		static constexpr int POISSON_RIGHT_MAX = 23;

	public:
    poisson(int size, ValueT* prob);
		~poisson();

		int left() const;
		int right() const;
		ValueT weight() const;
    ValueT lambda() const;

		ValueT operator()(int i) const;

    static int rightbound(ValueT lambda, ValueT eps);

    void set(ValueT lambda, int left, int right);

	private:
		ValueT m_lambda;
		int m_size;
		int m_left;
		int m_right;
		ValueT* m_prob;
		ValueT m_weight;

		static ValueT normalt(ValueT x);
		static ValueT normalq(ValueT p);
		ValueT set_prob();

	};

	template <typename ValueT>
	poisson<ValueT>::poisson(int size, ValueT* prob)
	: m_size(size), m_left(0), m_right(size-1), m_prob(prob) {}

	template <typename ValueT>
	poisson<ValueT>::~poisson() {}

  template <typename ValueT>
	inline int poisson<ValueT>::left() const {
		return m_left;
	}

	template <typename ValueT>
	inline int poisson<ValueT>::right() const {
		return m_right;
	}

	template <typename ValueT>
	inline ValueT poisson<ValueT>::weight() const {
		return m_weight;
	}

	template <typename ValueT>
	inline ValueT poisson<ValueT>::lambda() const {
		return m_lambda;
	}

	template <typename ValueT>
	inline ValueT poisson<ValueT>::operator()(int i) const {
		return m_prob[i - m_left];
	}

	template <typename ValueT>
	void poisson<ValueT>::set(ValueT lambda, int left, int right) {
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

	template <typename ValueT>
	ValueT poisson<ValueT>::normalt(ValueT x) {
		ValueT x2 = x*x;
		ValueT tmp = x;
		ValueT sum = 1 / tmp;
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

	template <typename ValueT>
	ValueT poisson<ValueT>::normalq(ValueT p) {
		const ValueT leps = log(p);
		assert(leps <= NORMALQ_UPPER_LOGP && leps >= NORMALQ_LOWER_LOGP);
		ValueT l = NORMALQ_LOWER_Q;
		ValueT u = NORMALQ_UPPER_Q;
		ValueT m = (l + u) / 2;
		ValueT fm = normalt(m) - leps;
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

	template <typename ValueT>
	int poisson<ValueT>::rightbound(ValueT lambda, ValueT eps) {
		int right;
		if (lambda == 0) {
			return 0;
		}
		if (lambda < POISSON_LAMBDA_MIN) {
			ValueT tmp = exp(-lambda);
			ValueT total = tmp;
			right = 0;
			for (int k=1; k<=POISSON_RIGHT_MAX; k++) {
				right++;
				tmp *= lambda / right;
				total += tmp;
				if (total + eps >= 1)
				break;
			}
		} else {
			ValueT z = normalq(eps);
			ValueT tmp = z + sqrt(4 * lambda - 1);
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

	template <typename ValueT>
	ValueT poisson<ValueT>::set_prob() {
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
		ValueT weight = 0;
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
