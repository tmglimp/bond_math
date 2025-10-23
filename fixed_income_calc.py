import pandas as pd
from math import pow
import config
from config import mvol
from datetime import datetime, timedelta
import logging
"""
SIA/FIA spot dirty formulas for:
 - No-arbitrage prices at market yield (BPrice, TPrice)
 - Accrued interest (aint)
 - Modified duration (MDur)
 - Macaulay duration (MacDur)
 - DV01 (PVBP)
 - Convexity (Cvx)
 - Basis overlay tails
 - Volcker Sensitivity
 - Volcker Stress
"""

def round_ytm(ytm):

    if pd.isnull(ytm):
        return None
    return round(ytm * 2) / 2.0

# ---------------- Date & Term Functions ----------------
def calculate_term(settlement_date_str, maturity_date_str, day_count_convention=365.25):

    settlement_date = datetime.strptime(settlement_date_str, '%Y%m%d')
    maturity_date = datetime.strptime(maturity_date_str, '%Y%m%d')
    days_to_maturity = (maturity_date - settlement_date).days
    term_in_years = days_to_maturity / day_count_convention
    return term_in_years

def compute_settlement_date(trade_date, t_plus=1):

    if isinstance(trade_date, str):
        trade_date = datetime.strptime(trade_date, '%Y%m%d')
    settlement_date = trade_date
    business_days_added = 0
    while business_days_added < t_plus:
        settlement_date += timedelta(days=1)
        if settlement_date.weekday() < 5:  # Monday=0, ..., Friday=4
            business_days_added += 1
    return settlement_date.strftime('%Y%m%d')

# ---------------- Yield and Price Functions ----------------
def accrual_period(begin, settle, next_coupon, day_count=1):

    if day_count == 1:
        L = datetime.strptime(str(begin), '%Y%m%d')
        S = datetime.strptime(str(settle), '%Y%m%d')
        N = datetime.strptime(str(next_coupon if next_coupon is not None else settle), '%Y%m%d')
        return (S - L).days / (N - L).days
    else:
        # 30/360 convention
        L = [int(begin[:4]), int(begin[4:6]), int(begin[6:8])]
        S = [int(settle[:4]), int(settle[4:6]), int(settle[6:8])]
        return (360 * (S[0] - L[0]) + 30 * (S[1] - L[1]) + S[2] - L[2]) / 180

def aint(cpn, period=2, begin=None, settle=None, next_coupon=None, day_count=1):

    v = accrual_period(begin, settle, next_coupon, day_count)
    return cpn / period * v

def BPrice(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    if term is None or yield_ is None:
        return None

    rounded_term = round_ytm(term)
    if rounded_term is None:
        return None

    T = int(rounded_term * period)  # Total coupon periods
    C = cpn / period
    Y = yield_ / period

    try:
        price = C * (1 - pow(1 + Y, -T)) / Y + 100 / pow(1 + Y, T)
    except ZeroDivisionError:
        price = None

    if cpn and begin and settle and next_coupon:
        ai = aint(cpn, period=2, begin=begin, settle=settle, next_coupon=next_coupon, day_count=1)
        price = price + ai

    return price

def TPrice(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    if term is None or yield_ is None:
        return None

    rounded_term = round_ytm(term)
    if rounded_term is None:
        return None
    T = int(rounded_term * period)  # Total coupon periods
    C = cpn / period
    Y = yield_ / period

    try:
        price = C * (1 - pow(1 + Y, -T)) / Y + 100 / pow(1 + Y, T)
    except ZeroDivisionError:
        price = None

    if begin and settle and next_coupon:
        ai = aint(cpn, period=2, begin=begin, settle=settle, next_coupon=next_coupon, day_count=day_count)
        price = price + ai
    return price/conv_factor

def MDur(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    if term is None or yield_ is None:
        return None

    rounded_term = round_ytm(term)
    if rounded_term is None:
        return None

    T = int(rounded_term * period)
    C = cpn / period
    Y = yield_ / period
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or P == 0:
        return None

    if begin and settle and next_coupon:
        v = accrual_period(begin, settle, next_coupon, day_count)
        P = pow(1 + Y, v) * P
        mdur = (-v * pow(1 + Y, v - 1) * C / Y * (1 - pow(1 + Y, -T))
                + pow(1 + Y, v) * (
                        C / pow(Y, 2) * (1 - pow(1 + Y, -T))
                        - T * C / (Y * pow(1 + Y, T + 1))
                        + (T - v) * 100 / pow(1 + Y, T + 1)))
    else:
        mdur = (C / pow(Y, 2) * (1 - pow(1 + Y, -T))) + (T * (100 - C / Y) / pow(1 + Y, T + 1))
    return mdur / (period * P)

def MacDur(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if mdur is None:
        return None
    return mdur * (1 + yield_ / period)

def DV01(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or mdur is None:
        return None
    return round((mdur) * P * 0.001, 6)

def fut_DV01(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or mdur is None:
        return None
    return round((mdur) * P * 0.001, 6)/ conv_factor

def DV01minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ - .0001
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or mdur is None:
        return None
    return round((mdur) * P * 0.001, 6)

def fut_DV01minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ - .0001
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or mdur is None:
        return None
    return round((mdur) * P * 0.001, 6)/ conv_factor

def DV01plus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or mdur is None:
        return None
    return round((mdur + cvx) * P * 0.001, 6)

def DV10(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_+.001
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV10(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_+.001
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def DV10minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_-.001
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV10minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_-.001
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6) /conv_factor

def DV50(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ + .005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV50(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ + .005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/ conv_factor

def DV50minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ - .005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV50minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ - .005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6) / conv_factor

def DV100(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ + .01
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV100(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ + .01
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/ conv_factor

def DV100minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ - .01
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_DV100minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ - .01
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6) / conv_factor

def Cvx(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    if term is None or yield_ is None:
        return None
    rounded_term = round_ytm(term)
    if rounded_term is None:
        return None
    T = int(rounded_term * period)
    C = cpn / period
    Y = yield_ / period
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or P == 0:
        return None
    v = accrual_period(begin, settle, next_coupon, day_count) if (begin and settle and next_coupon) else 0
    dcv = (
            -v * (v - 1) * pow(1 + Y, v - 2) * C / Y * (1 - pow(1 + Y, -T))
            - 2 * v * pow(1 + Y, v - 1) * (C / pow(Y, 2) * (1 - pow(1 + Y, -T)) - T * C / (Y * pow(1 + Y, T + 1)))
            - pow(1 + Y, v) * (
                    -C / pow(Y, 3) * (1 - pow(1 + Y, -T)) +
                    2 * T * C / (pow(Y, 2) * pow(1 + Y, T + 1)) +
                    T * (T + 1) * C / (Y * pow(1 + Y, T + 2))
            )
            + (T - v) * (T + 1) * 100 / pow(1 + Y, T + 2 - v)
    )
    return dcv / (P * period ** 2)

def fut_Cvx(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    if term is None or yield_ is None:
        return None
    rounded_term = round_ytm(term)
    if rounded_term is None:
        return None
    T = int(rounded_term * period)
    C = cpn / period
    Y = yield_ / period
    P = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if P is None or P == 0:
        return None
    v = accrual_period(begin, settle, next_coupon, day_count) if (begin and settle and next_coupon) else 0
    dcv = (
            -v * (v - 1) * pow(1 + Y, v - 2) * C / Y * (1 - pow(1 + Y, -T))
            - 2 * v * pow(1 + Y, v - 1) * (C / pow(Y, 2) * (1 - pow(1 + Y, -T)) - T * C / (Y * pow(1 + Y, T + 1)))
            - pow(1 + Y, v) * (
                    -C / pow(Y, 3) * (1 - pow(1 + Y, -T)) +
                    2 * T * C / (pow(Y, 2) * pow(1 + Y, T + 1)) +
                    T * (T + 1) * C / (Y * pow(1 + Y, T + 2))
            )
            + (T - v) * (T + 1) * 100 / pow(1 + Y, T + 2 - v)
    )
    return dcv / (P * period ** 2)/conv_factor

def sensitivity22(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ + .0002
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def sensitivity22minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ - .0002
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def sensitivity55(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ + .0005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def sensitivity55minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1):
    yield_ = yield_ - .0005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_sensitivity22(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    yield_ = yield_ + .0002
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def fut_sensitivity22minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    yield_ = yield_ - .0002
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def fut_sensitivity55(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,conv_factor=None):
    yield_ = yield_ + .0005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def fut_sensitivity55minus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, conv_factor=None):
    yield_ = yield_ - .0005
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def sensitivityMKT(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, mvol = None):
    yield_ = yield_ + mvol
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def sensitivityMKTminus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, mvol = None):
    yield_ = yield_ - mvol
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)

def fut_sensitivityMKT(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, mvol = None,conv_factor=None):
    yield_ = yield_ + mvol
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def fut_sensitivityMKTminus(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1, mvol = None,conv_factor=None):
    yield_ = yield_ - mvol
    mdur = MDur(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    #cvx = Cvx(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    if yield_ is None or mdur is None:
        return None
    return round((mdur) * 0.001, 6)/conv_factor

def appx_duration(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,delta_y=0.0001):
    if yield_ is None:
        return None
    price = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    p_bp_up = BPrice(cpn, term, (yield_ + delta_y), period, begin, settle, next_coupon, day_count)
    p_bp_dn = BPrice(cpn, term, (yield_ - delta_y), period, begin, settle, next_coupon, day_count)
    if price is None or price == 0:
        return None
    return (p_bp_dn - p_bp_up) / (2 * price * delta_y)

def appx_convexity(cpn, term, yield_, period=2, begin=None, settle=None, next_coupon=None, day_count=1,delta_y=0.0001):
    if yield_ is None:
        return None
    price = BPrice(cpn, term, yield_, period, begin, settle, next_coupon, day_count)
    p_bp_up = BPrice(cpn, term, (yield_ + delta_y), period, begin, settle, next_coupon, day_count)
    p_bp_dn = BPrice(cpn, term, (yield_ - delta_y), period, begin, settle, next_coupon, day_count)
    if price is None or price == 0:
        return None
    return (p_bp_up + p_bp_dn - 2 * price) / (price * delta_y ** 2)

def sia_implied_repo(fut_price, dirty_price, cf, days):
    adj_fut = fut_price * cf
    return (((adj_fut - dirty_price) / dirty_price)-1) * (365 / days)

def sia_gross_basis(fut_price, cf, dirty_price):
    return fut_price * cf - dirty_price

def sia_convexity_yield(dirty_price, coupon, days):
    return ((coupon / dirty_price) * (days / 365)) * (365 / days)

def sia_carry(gross_basis, implied_repo, dirty_price, days):
    financing_cost = dirty_price * implied_repo * days / 365
    return gross_basis - financing_cost

def sia_net_basis(gross_basis, carry):
    return gross_basis + carry

def fut_tail(A_FUT_DV01, A_MULT, B_FUT_DV01, B_MULT):
    return (-1 *
        (((A_FUT_DV01 * A_MULT) - (B_FUT_DV01 * B_MULT)) / (B_FUT_DV01 * B_MULT))
        if (A_FUT_DV01 * A_MULT) > (B_FUT_DV01 * B_MULT)
        else (((B_FUT_DV01 * B_MULT) - (A_FUT_DV01 * A_MULT)) / (A_FUT_DV01 * A_MULT))
    )

def fwd_fut_tail(A_FUT_DV01, A_FWD_DV01, A_MULT, B_FUT_DV01, B_FWD_DV01, B_MULT):
    return (-1 *
        ((((A_FUT_DV01+A_FWD_DV01) * A_MULT) - ((B_FUT_DV01+B_FWD_DV01) * B_MULT)) /
                ((B_FUT_DV01+B_FWD_DV01) * B_MULT))
        if (A_FUT_DV01 * A_MULT) > (B_FUT_DV01 * B_MULT)
        else ((((B_FUT_DV01+B_FWD_DV01) * B_MULT) - ((A_FUT_DV01+A_FWD_DV01) * A_MULT)) /
                ((A_FUT_DV01+A_FWD_DV01) * A_MULT))
    )

def calculate_ytm(market_price, face_value, coupon_rate, time_to_maturity, periods_per_year=2, n_digits=4):
    # Convert market_price from percent to actual price in currency terms.
    market_price = market_price / 100.0 * face_value
    coupon_rate = coupon_rate / 100.0
    coupon_payment = face_value * coupon_rate / periods_per_year

    def bond_price(ytm):
        pv = 0
        T = int(time_to_maturity * periods_per_year)
        for t in range(1, T + 1):
            pv += coupon_payment / (1 + ytm / periods_per_year) ** t
        pv += face_value / (1 + ytm / periods_per_year) ** T
        return pv

    ytm_guess = coupon_rate
    tolerance = 1e-6
    max_iterations = 1000
    ytm = ytm_guess

    for _ in range(max_iterations):
        price_at_ytm = bond_price(ytm)
        delta_ytm = 1e-5
        p_bp_up = bond_price(ytm + delta_ytm)
        p_bp_dn = bond_price(ytm - delta_ytm)
        price_derivative = (p_bp_up - p_bp_dn) / (2 * delta_ytm)

        if abs(price_derivative) < 1e-8:
            return round(ytm, n_digits)

        ytm_new = ytm - (price_at_ytm - market_price) / price_derivative
        if abs(ytm_new - ytm) < tolerance:
            return round(ytm_new, n_digits)
        ytm = ytm_new
    print("YTM calculation did not converge within the maximum number of iterations.")
    return ytm

def _log_kpis(prefix: str = ""):
    logging.info(
        f"{prefix} KPIs :: "
        f"FUTFV:{getattr(config,'FUT_TPRICE', None)} | "
        f"CTDFV:{getattr(config,'CTD_BPRICE', None)} | "
        f"CTDMDUR:{getattr(config,'CTD_MDUR', None)} | "
        f"CTDMACDUR:{getattr(config,'CTD_MACDUR', None)} | "
        f"CTDCVX:{getattr(config,'CTD_CVX', None)} | "
        f"CTDDV01:{getattr(config,'CTD_DV01', None)} | "
        f"FUTCVX:{getattr(config,'FUT_CVX', None)} | "
        f"FUTDV01:{getattr(config,'FUT_DV01', None)} | "
        f"FUTDV01-:{getattr(config,'FUT_DV01_MINUS', None)} | "
        f"FUTDV10:{getattr(config,'FUT_DV10', None)} | "
        f"FUTDV10-:{getattr(config,'FUT_DV10_MINUS', None)} | "
        f"FUTDV50:{getattr(config,'FUT_DV50', None)} | "
        f"FUTDV50-:{getattr(config,'FUT_DV50_MINUS', None)} | "
        f"FUTDV100:{getattr(config,'FUT_DV100', None)} | "
        f"FUTDV100-:{getattr(config,'FUT_DV100_MINUS', None)} | "
        f"FUTDV22:{getattr(config,'FUT_DV22', None)} | "
        f"FUTDV22-:{getattr(config,'FUT_DV22_MINUS', None)} | "
        f"FUTDV55:{getattr(config,'FUT_DV55', None)} | "
        f"FUTDV55-:{getattr(config,'FUT_DV55_MINUS', None)}"
    )