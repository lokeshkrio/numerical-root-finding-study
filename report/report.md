# **1. Introduction**

Finding the roots of nonlinear equations is an important problem in numerical mathematics. Many equations cannot be solved easily by using analytical methods, so numerical algorithms are used to obtain approximate solutions to the equations.

This project compares four\* numerical root-finding algorithms: **Bisection, Regula Falsi, Newton-Raphson, Secant**. The algorithms are implemented in Python and tested on different mathematical functions.

The performance of each algorithm is compared on the basis of the number of iterations, execution time, residual error, and accuracy of the calculated root.

The main objective is to understand the efficiency, accuracy, and limitations of these algorithms through practical experiments.

\*will add Müller in future

# **2. Algorithms**


I. Bracketing Methods

*Require an initial interval $[a,b]$ where a sign change occurs ($f(a)\cdot f(b)<0$). They guarantee convergence if $f(x)$ is continuous.*

**1. Bisection Method**

* **Concept:** Continually halves the interval containing the root.
* **Formula:** Midpoint $c=\frac{a+b}{2}$.
* **Sign Check:** If $f(a)\cdot f(c)<0$, root lies in $[a,c]$ $\rightarrow$ set $b=c$. Else, set $a=c$.
* **Convergence Rate:** **Linear** (slow, but 100% reliable).
* **Recurrence Relation:** The next estimate ${c}_{k}$ is computed as the arithmetic mean:
  			${c}_{k}=\frac{{a}_{k}+{b}_{k}}{2}$
* **Pro/Con:** Requires no derivative; cannot directly locate even roots (where the function touches the x-axis without crossing).

**2. Regula Falsi (False Position Method)**

* **Concept:** Instead of using the midpoint, it connects $(a,f(a))$ and $(b,f(b))$ with a straight line. The x-intercept of this line serves as the next estimate.
* **Formula:** $c=\frac{a\cdot f(b)-b\cdot f(a)}{f(b)-f(a)}$ *(or equivalently $c=b-\frac{f(b)(b-a)}{f(b)-f(a)}$)*.
* **Recurrence Relation:** The estimate ${c}_{k}$ is the zero-crossing of the secant line intersecting $({a}_{k},f({a}_{k}))$ and $({b}_{k},f({b}_{k}))$:
  ${c}_{k}={b}_{k}-\frac{f({b}_{k})({b}_{k}-{a}_{k})}{f({b}_{k})-f({a}_{k})}=\frac{{a}_{k}f({b}_{k})-{b}_{k}f({a}_{k})}{f({b}_{k})-f({a}_{k})}$
* **Convergence Rate:** **Linear** (typically faster than bisection, but can become extremely slow if one endpoint stays fixed).
* **Pro/Con:** Better speed than bisection for most well-behaved functions, but prone to one-sided stagnation.

---

II. Open Methods

*Use local information (points or derivatives) near a single guess. They do not require the root to be bracketed, but they can diverge if initial guesses are poor.*

**1. Newton-Raphson Method**

* **Concept:** Uses the tangent line of the function at the current point to project the next estimate onto the x-axis.
* **Formula:** ${x}_{n+1}={x}_{n}-\frac{f({x}_{n})}{f'({x}_{n})}$.
* **Convergence Rate:** **Quadratic** (very fast near a simple root).
* **Pro/Con:** Fast execution; fails if $f'({x}_{n})=0$ (horizontal tangent) or if the operation cycles infinitely.

**2. Secant Method**

* **Concept:** Eliminates Newton’s requirement for a derivative by approximating $f'(x)$ using a secant line through two historic points.
* **Formula:** ${x}_{n+1}={x}_{n}-f({x}_{n})\cdot \frac{{x}_{n}-{x}_{n-1}}{f({x}_{n})-f({x}_{n-1})}$
* **Convergence Rate:** **Superlinear** (Order $\approx 1.618$, the golden ratio).
* **Pro/Con:** No analytical derivative needed; requires two initial points to start, and can still diverge.

**3. Müller’s Method**

* **Concept:** Generalisation of the Secant method. Instead of drawing a straight line through two points, it fits a **quadratic parabola** through three points $({x}_{0},{x}_{1},{x}_{2})$.
* **Formula:** Solves for the root of the fitted polynomial $a(x-{x}_{2}{)}^{2}+b(x-{x}_{2})+c=0$:
  ${x}_{3}={x}_{2}-\frac{2c}{b\pm \sqrt{{b}^{2}-4ac}}$

  *(Note: Choose the $\pm$ sign to match the sign of $b$ to maximise the denominator and ensure numerical stability)*.
* **Convergence Rate:** **Superlinear** (Order $\approx 1.84$).
* **Pro/Con:** Can discover **complex/imaginary roots** even with real initial guesses (due to the square root); algebra is more complex to implement.

---

III. Hybrid Methods

**Brent’s Method**

* **Concept:** Combines the rock-solid reliability of bracketing with the high speed of open methods. It tracks three bounds and swaps mechanisms dynamically based on internal criteria.
* **Components Used:**
  * **Bisection** (Fallback safety net)
  * **Secant Method / Linear Interpolation** (Used when tracking 2 points)
  * **Inverse Quadratic Interpolation (IQI)** (Fits $x$ as a quadratic function of $y$ using 3 points; computes next root by evaluating at $y=0$).
* **How it Works:**
  * Attempts **IQI** or **Secant** steps first to achieve fast, superlinear convergence.
  * Runs strict safety checks at every step: If the fast interpolation attempts to step out of bounds or isn't shrinking the interval quickly enough, it immediately **forces a Bisection step**.
* **Convergence Rate:** Behaves **superlinearly** on smooth functions, but defaults to **linear** in the worst-case scenario.
* **Pro/Con:** Globally convergent like bisection, yet fast like open methods. The industry standard for single-variable black-box root finding

## **2.5 Summary**

| Algorithm | Initial Requirements | Derivative | Approach | Key Characteristic |
| ----- | ----- | ----- | ----- | ----- |
| **Bisection** | Two endpoints | Not required | Bracketing | Reliable, but relatively slow |
| **Regula Falsi** | Two endpoints | Not required | Bracketing | Uses linear interpolation |
| **Newton-Raphson** | One initial guess | Required | Open | Very fast near a simple root |
| **Secant** | Two initial guesses | Not required | Open | Approximates the derivative |
| **Müller** | Three initial points | Not required | Open | Uses quadratic interpolation |
| **Brent** | Two endpoints forming a valid bracket | Not required | Hybrid | Combines robustness with fast convergence |

*Note:* Bracketing methods maintain an interval containing a root. Open methods generate successive estimates without maintaining such an interval. Brent’s method combines bracketing, interpolation, and bisection safeguards.

# **3. Implementation**

The algorithms were implemented in Python using separate functions for each method. Experiments were conducted on different mathematical functions to compare their performance.

## **3.1 Experimental Setup**

* **Language:** Python
* **Algorithms:** Bisection, Regula Falsi, Newton–Raphson, Secant
* **Test functions:** 11 mathematical functions
* **Tolerance:** (10^{-6})
* **Maximum iterations:** 1000
* **Repeated trials:** 20 per algorithm–function combination

## **3.2 Performance Metrics**

The algorithms were evaluated using the following metrics:

* **Number of iterations:** Measures how many iterations an algorithm requires.
* **Execution time:** Measures the time taken to calculate the root.
* **Residual error:** Measures how close the calculated root is to satisfying (f(x)=0).
* **Root error:** Measures the difference between the calculated root and the known reference root.

## **3.3 Validation and Error Handling**

The implementation includes checks for invalid inputs, division by zero, non-finite values, and iteration limits. The results are saved in a CSV file for analysis and visualization.

Python's Matplotlib and Pandas libraries are used to analyze the results and generate comparison plots.

# **4. Experimental Results and Analysis**

The four algorithms were tested on 11 mathematical functions with a tolerance of $1{0}^{-6}$, a maximum of 1000 iterations, and 20 repeated trials.

## **4.1 Performance Comparison**

| Algorithm | Average iterations | Functions matching target root |
| ----- | ----- | ----- |
| Bisection | 20.18 | 11/11 |
| Regula Falsi | 13.36 | 11/11 |
| Newton–Raphson | 3.82 | 11/11 |
| Secant | 8.00 | 10/11 |

## **4.2 Observations**

* **Newton-Raphson** required fewest iterations on average.
* **Regula Falsi** was generally more efficient than Bisection in terms of iteration count.
* **Bisection** consistently found the target roots but sometimes had residuals above the specified tolerance.
* **Secant** performed well on most functions but converged to a different root for one test function.
* Execution times were extremely small, so iteration counts and accuracy provide more useful comparisons than tiny timing differences.


## **4.3 Conclusion from Results**

Newton-Raphson was the most efficient method in terms of iteration count.

Bisection offered a reliable bracketing approach, while Regula Falsi was an alternative bracketing approach.

The Secant method performed well overall but showed that finding a root does not necessarily mean finding the particular target root.

These results highlight the importance of considering accuracy, convergence behaviour, and the properties of each function when selecting a numerical root-finding method.

# **5. Conclusion and Future Improvements**

## **5.1 Conclusion**

This project compared four numerical root-finding algorithms using Python. The experiments showed that Newton-Raphson required the fewest iterations on average, while Bisection provided reliable bracketing. Regula Falsi offered an alternative bracketing method, and the Secant method performed well on most test functions.

The study demonstrates that no single algorithm is ideal for every function. The choice of method depends on convergence behaviour, accuracy requirements, and the availability of derivatives or a suitable initial interval.

## **5.2 Future Improvements**

The project can be extended by:

* Testing more functions, including repeated roots.
* Studying the effect of different initial guesses and tolerances.
* Adding other methods, such as Müller's method.
* Comparing convergence rates across different classes of functions.
