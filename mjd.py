import math
from scipy.stats import norm

JF_INTENSITY = 0.25
JF_DRIFT = -0.015
JF_VOLATILITY = 0.04

VOLATILITY = 0.15
PRICE = 50
STRIKE_PRICE = 50

RF_RATE = 0.02
DIV = 0
TAO = 0.25

def get_d_plus(x, K, r, q, sigma, tao):
    num = math.log(x/K) + (r-q+((sigma**2)/2)) * tao
    denom = sigma * tao**0.5
    return num / denom

def get_N_d_plus(x, K, r, q, sigma, tao):
    d_plus = get_d_plus(x, K, r, q, sigma, tao)
    return norm.cdf(d_plus, loc=0, scale=1)

def get_d_minus(x, K, r, q, sigma, tao):
    plus = get_d_plus(x, K, r, q, sigma, tao)
    return plus - (sigma * tao**0.5)

def get_N_d_minus(x, K, r, q, sigma, tao):
    d_minus = get_d_minus(x, K, r, q, sigma, tao)
    return norm.cdf(d_minus, loc=0, scale=1)

def get_bsm_call_price(price, K, sigma, r, tao, q):
    first = price * math.exp(-1 * q * tao) * get_N_d_plus(price, K, r, q, sigma, tao)
    second = K * math.exp(-1 * r * tao) * get_N_d_minus(price, K, r, q, sigma, tao)
    return first - second

def get_mjd_call_price(intensity, tao, price, K, r, q, sigma, sigma_j, mu_j, N=50):
    ret = 0
    a = math.exp(mu_j + ((sigma_j**2)/2)) - 1
    outside_factor = math.exp(-1 * intensity * (1+a) * tao)
    for i in range(N):
        s_i = (sigma**2) + ((i/tao) * ((sigma_j**2)/2))
        r_i = r - (intensity * a) + ((i/tao)*(math.log(1+a)))
        bsm_i = get_bsm_call_price(price, K, s_i**0.5, r_i, tao, q)
        other_i = ((intensity * (1+a) * tao)**i) / (math.factorial(i))
        ret = ret + (bsm_i * other_i * outside_factor)
    return ret

price = get_mjd_call_price(JF_INTENSITY, TAO, PRICE, STRIKE_PRICE, RF_RATE, DIV,
                      VOLATILITY, JF_VOLATILITY, JF_DRIFT)
print(price)
