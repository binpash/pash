// a suite groups assertions such that failure halts only the suite.

const suite = async (msg, f) => {
    try {
        return await Promise.resolve(f());
    } catch (e) {
        if (e.suite) {
            delete e.suite;
            console.log(msg, '\n', e);
        }
        else throw e;
    }
};

const assert = (msg, v) => {
    if (!v) {
        const e = new Error(msg);
        e.suite = true;
        throw e;
    }
};

const eqv = (a, b, e = 0.001) => Math.abs(a - b) < e;

const assertEqv = (msg, a, b, e = 0.001) =>
      assert((msg && `${msg}: ` + `${a} ≈ ${b} (within ${e})`), eqv(a, b, e));

module.exports = {
    assert,
    suite,
    eqv,
    assertEqv,
};
