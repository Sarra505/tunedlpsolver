#pragma once

#include "SolverAlgorithm.hpp"

struct Random
{
    /// The state
    uint64_t randomState;

    /// Constructor
    explicit Random(uint64_t seed = 0) : randomState(1152921504606847009ull ^ (seed * 1099511627791ull)) {}

    /// Add constant
    static constexpr uint64_t addConstant = 0xa0761d6478bd642full;
    /// Xor constant
    static constexpr uint64_t xorConstant = 0xe7037ed1a0b428dbull;

    /// Get the next value
    uint64_t next()
    {
        randomState += addConstant;
        // Multiply fold
        auto m = static_cast<unsigned __int128>(randomState) * static_cast<unsigned __int128>(randomState ^ xorConstant);
        return static_cast<uint64_t>(m >> 64) ^ static_cast<uint64_t>(m);
    }
    /// Get the next value in range [0, s)
    /// Based on Lemire, Fast Random Integer Generation in an Interval https://arxiv.org/abs/1805.10941
    uint64_t nextRange(uint64_t s) noexcept
    {
        uint64_t x = next();
        auto m = static_cast<unsigned __int128>(x) * s;
        auto l = static_cast<uint64_t>(m);
        if (l < s) [[unlikely]]
        { // unlikely as we assume s to be smaller than UINT64_MAX
            uint64_t t = -s % s;
            while (l < t)
            {
                x = next();
                m = static_cast<unsigned __int128>(x) * s;
                l = static_cast<uint64_t>(m);
            }
        }
        return m >> 64;
    }
    /// Get the next value in range [start, end)
    uint64_t nextRange(uint64_t start, uint64_t end)
    {
        return start + nextRange(end - start);
    }
    /// Result type, required for C++ Concept UniformRandomBitGenerator
    using result_type = uint64_t;
    /// Call. Alias for next(), required for C++ Concept UniformRandomBitGenerator
    auto operator()() { return next(); }
    /// Smallest value returned by Random, required for C++ Concept UniformRandomBitGenerator
    static constexpr uint64_t min() { return 1; }
    /// Largest value returned by Random, required for C++ Concept UniformRandomBitGenerator
    static constexpr uint64_t max() { return ~0; }
};
//---------------------------------------------------------------------------
// Approximate packing linear program solver with (1 / eps) multiplicative error
// Optimizes max_{x >= 0}(sum(x) : Ax <= b)
class RevisedSimplexMPFI
{
public:
    /// A buffer of aligned floats
    struct alignas(64) FloatBuffer
    {
        /// The number of floats in a cache line
        static constexpr size_t Size = 64 / sizeof(float);
        /// The floats
        float buffer[Size]{};
    };

    /// The sparse coefficient matrix M
    std::vector<Entry> matrix;
    /// Compact sparse M^T coefficients
    std::vector<float> coefs;
    /// step number for last solved LP
    unsigned lastStepCount;
    /// Compact sparse M^T indices
    std::vector<unsigned> rules;
    /// Compact sparse M^T index offsets per column
    std::vector<unsigned> colOffsets;
    /// Storage for solving
    std::vector<FloatBuffer> storage;
    /// Additional buffer for computations
    std::vector<unsigned> buffer;
    /// MPFU coefficients
    std::vector<FloatBuffer> mpfu_a;
    /// MPFU leaving variable indices
    std::vector<unsigned> mpfu_p;
    /// Basic and nonbasic column indices
    std::vector<unsigned> basic, nonBasic;
    /// The number of rules
    unsigned m = 0;

    auto getCoefs(unsigned i) { return coefs.data() + colOffsets[i]; };
    auto getRules(unsigned i) { return rules.data() + colOffsets[i]; };
    auto getSz(unsigned i) { return colOffsets[i + 1] - colOffsets[i]; }

    /// Prepare for computation
    void prepare(unsigned n);
    /// Solve the given LP using simplex
    float solveSimplex(unsigned n, unsigned stepLimit = ~0u);
};
//---------------------------------------------------------------------------
void RevisedSimplexMPFI::prepare(unsigned n)
// Prepare for computation
{
    // Build optimized column representation for fast operations
    colOffsets.clear();
    colOffsets.resize(n + 1);
    rules.resize(matrix.size());
    coefs.resize(matrix.size());

    sort(matrix.begin(), matrix.end(), [](auto &a, auto &b)
         { return tie(a.variable, a.rule) < tie(b.variable, b.rule); });

    for (size_t i = 0; i < matrix.size(); i++)
    {
        assert(matrix[i].rule < m);
        assert(matrix[i].variable < n);
        rules[i] = matrix[i].rule;
        coefs[i] = matrix[i].coef;
        colOffsets[matrix[i].variable + 1]++;
    }

    for (unsigned i = 1; i <= n; i++)
        colOffsets[i] += colOffsets[i - 1];
}
//---------------------------------------------------------------------------
float RevisedSimplexMPFI::solveSimplex(unsigned n, unsigned stepLimit)
// Solve the given LP using simplex
{
    if (n == 0)
        return 0.0;
    if (m == 0)
        return numeric_limits<float>::infinity();

    prepare(n);

    // We use the revised simplex algorithm (according to Chvatal, "Linear Programming") since it is fast
    // with MPFU (Middle Product Form Update, due to Huangfu & Hall, "Novel update techniques for the revised simplex method") since it is a good balance between fast and simple to implement.

    // Our basis matrix is initially the identity matrix
    // It's PLU factorization is of course also trivially the identity matrix
    // Since we assume our simplex algorithm is only executed for a small number of steps,
    //   we never refactorize, and keep the identity assumption throughout
    //   this makes our implementation very simple (and fast for a small number of steps) as many formulas become simpler
    // For example, as we assume U = I, the MPFU variable e^~_p := u_p = e_p = (0, 0, ..., 0, 1, 0, 0, ..., 0), where only the index p is set to 1

    // The MPFU steps need two variables
    // The vector a = (a^~_p - u_p)/mu where u_p = e_p and mu = 1 + e_p @ (a^~_p - u_p)
    // And the unsigned p, indicating the leaving variable at that step
    size_t mpfu_aSize = 0;

    mpfu_a.clear();
    mpfu_p.clear();
    // We do not allow the vector to be empty to avoid undefined behaviour
    mpfu_a.emplace_back();

    basic.clear();
    nonBasic.clear();
    basic.reserve(m);
    nonBasic.reserve(n);
    for (unsigned i = 0; i < n; i++)
        nonBasic.push_back(i);
    for (unsigned j = n; j < n + m; j++)
        basic.push_back(j);

    size_t nCount = 0;
    size_t mCount = 2;
    storage.resize((n + FloatBuffer::Size - 1) / FloatBuffer::Size * nCount + (m + FloatBuffer::Size - 1) / FloatBuffer::Size * mCount);
    unsigned lastAlloc = 0;
    auto alloc = [&](unsigned sz)
    {
        auto *ptr = storage[lastAlloc].buffer;
        lastAlloc += (sz + FloatBuffer::Size - 1) / FloatBuffer::Size;
        assert(lastAlloc <= storage.size());
        return ptr;
    };
    float *__restrict xb = alloc(m);
    float *__restrict y = alloc(m);

    for (unsigned j = 0; j < m; j++)
        xb[j] = 1.0f;

    auto applyMPFUForwardStep = [&](float *__restrict a, float *__restrict b, unsigned p)
    {
        float y = b[p];
        for (unsigned j = 0; j < m; j++)
            b[j] -= y * a[j];
    };
    auto applyMPFUBackwardStep = [&](float *__restrict a, float *__restrict b, unsigned p)
    {
        float y = 0.0f;
        for (unsigned j = 0; j < m; j++)
            y += b[j] * a[j];
        b[p] -= y;
    };
    // Solve B @ x = b for x
    auto applyMPFUForward = [&](float *__restrict b)
    {
        float *st = mpfu_a[0].buffer;
        float *en = mpfu_a[0].buffer + mpfu_aSize;
        unsigned *p = mpfu_p.data();
        for (auto *it = st; it != en; it += m, ++p)
            applyMPFUForwardStep(it, b, *p);
    };
    // Solve x @ B = B^T @ x = b for x
    auto applyMPFUBackward = [&](float *__restrict b)
    {
        float *en = mpfu_a[0].buffer - m;
        float *st = mpfu_a[0].buffer + mpfu_aSize - m;
        unsigned *p = mpfu_p.data() + mpfu_p.size() - 1;
        for (auto *it = st; it != en; it -= m, --p)
            applyMPFUBackwardStep(it, b, *p);
    };

    Random rng;

    unsigned step = 0;
    bool foundSolution = false;
    for (; step < stepLimit; step++)
    {
        // We want to find a variable with a positive coefficient
        // Increasing this variable should increase the objective
        // y := the coefficients for basic variables
        for (unsigned ind = 0; ind < m; ind++)
            y[ind] = basic[ind] < n ? 1.0f : 0.0f;
        applyMPFUBackward(y);

        // Should we pick the first valid index as entering to prevent cycling?
        bool pickFirst = rng.next() & 0b11;

        unsigned enteringInd = ~0u;
        float bestDot = 1e-5f;
        // Find the component of c_n - y @ N that is positive
        // If no positive column exists, then we found the optimal
        // We choose the column with the highest coefficient at the risk of cycling
        for (unsigned i = 0; i < n; i++)
        {
            unsigned col = nonBasic[i];
            float dot;
            if (col < n)
            {
                const float *__restrict cs = getCoefs(col);
                auto rs = getRules(col);
                auto sz = getSz(col);

                dot = 1.0f;
                for (unsigned k = 0; k < sz; k++)
                    dot -= y[rs[k]] * cs[k];
            }
            else
            {
                dot = -y[col - n];
            }
            if (dot > bestDot)
            {
                bestDot = dot;
                enteringInd = i;
                if (pickFirst)
                    break;
            }
        }

        // Are we at the optimum?
        if (enteringInd == ~0u)
        {
            foundSolution = true;
            break;
        }

        // y := the entering column
        for (unsigned j = 0; j < m; j++)
            y[j] = 0.0f;
        {
            unsigned col = nonBasic[enteringInd];
            if (col < n)
            {
                const float *__restrict cs = getCoefs(col);
                auto rs = getRules(col);
                auto sz = getSz(col);

                for (unsigned k = 0; k < sz; k++)
                    y[rs[k]] = cs[k];
            }
            else
            {
                y[col - n] = 1.0f;
            }
        }
        applyMPFUForward(y);

        // Find the largest value of t such that xb >= t * y
        // If no such t exists, the problem is unbounded
        float t = numeric_limits<float>::infinity();
        unsigned leavingInd = ~0u;
        for (unsigned j = 0; j < m; j++)
        {
            assert(xb[j] >= 0.0f);
            if (y[j] <= 0.0f)
                continue;
            float div = xb[j] / y[j];
            if (div < t)
            {
                leavingInd = j;
                t = div;
            }
        }

        // Check if problem is unbounded
        if (leavingInd == ~0u)
            return numeric_limits<float>::infinity();

        for (unsigned j = 0; j < m; j++)
        {
            // Prevent xb from turning negative due to slight numerical errors
            xb[j] = max(xb[j] - y[j] * t, 0.0f);
        }

        // Compute the vector a = (a^~_p - u_p)/mu where u_p = e_p and mu = 1 + e_p @ (a^~_p - u_p) and p = leavingInd
        float mu = y[leavingInd];
        float invMu = 1.0f / mu;
        y[leavingInd] -= 1.0f;
        for (unsigned j = 0; j < m; j++)
            y[j] *= invMu;

        if (mpfu_a.size() * FloatBuffer::Size < mpfu_aSize + m)
            mpfu_a.resize((mpfu_aSize + m + FloatBuffer::Size - 1) / FloatBuffer::Size);

        {
            float *__restrict base = mpfu_a[0].buffer + mpfu_aSize;
            for (unsigned j = 0; j < m; j++)
                base[j] = y[j];
            mpfu_aSize += m;
            mpfu_p.push_back(leavingInd);
        }

        swap(basic[leavingInd], nonBasic[enteringInd]);
        xb[leavingInd] = t;
    }

    if (!foundSolution)
        return -1.0f;

    float result = 0.0f;
    for (unsigned ind = 0; ind < m; ind++)
        result += basic[ind] < n ? xb[ind] : 0.0f;

    lastStepCount = step + 1;
    return result;
}
