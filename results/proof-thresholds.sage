import functools

@functools.cache
def calculate_V_list(m,t):
    V_list = [1]
    for d in range(1, t*m+1):
        V = V_list[-1] + binomial(m*t,d)
        V_list.append(V)
    return V_list

@functools.cache
def success_probability_from_error(m, t, error):
    V_list = calculate_V_list(m, t)
    S_1 = 0
    S_2 = 0
    for d in range(t*m+1):
        q_1 = RR( binomial(m*t,d) * error^d * (1-error)^(m*t-d) )
        q_2 = RR( (1 - V_list[d] / 2^(t*m))^(2^m) )               
        S_1 += (q_2)^(2^m) * q_1
        S_2 += (q_2) * q_1
    return S_1^(t+1) * S_2^(m*t-(t+1))
    
def proof_threshold(m, t, error):
    p1 = success_probability_from_error(m, t, error)
    p2 = success_probability_from_error(m, t, error + 0.001)
    try:
        assert(p1 >= 0.5)
        assert(p2 < 0.5)
    except AssertionError as exc:
        extra = 3
        es = [error + (x - extra) * 0.001 for x in range(2*extra)]
        for e in es:
            p = success_probability_from_error(m, t, e)
            print("{};{};{};{}".format(m, t, e, p))
        raise exc
    print("Threshold error probabilty for m={}, t={}: {}".format(m, t, error))

proof_threshold(9,20,0.261)
proof_threshold(10,50,0.340)
proof_threshold(12,64,0.360)
proof_threshold(13,96,0.384)
proof_threshold(13,119,0.395)
proof_threshold(13,128,0.398)
