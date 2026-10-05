import numpy as np

from fdm_solver import run_fdm


# --------------------------------------------------------------------------------------
# Common physical parameters
# --------------------------------------------------------------------------------------

V_UP = 5.0e-6
V_DOWN = 5.0e-6

L = 0.05
D = 0.0254

POROSITY = 0.10
K = 1.0e-17
MU = 1.8e-5
B = 1.0e5

N = 40


# --------------------------------------------------------------------------------------
# Test 1: Uniform pressure
# --------------------------------------------------------------------------------------

def test_uniform_pressure():
    """
    If upstream and downstream pressures are initially identical,
    no pressure gradient exists and the entire system should remain
    unchanged to numerical roundoff.
    """

    P_initial = 1.0e6

    result = run_fdm(
        V_up=V_UP,
        V_down=V_DOWN,
        L=L,
        D=D,
        porosity=POROSITY,
        k=K,
        mu=MU,
        b=B,
        P_up_initial=P_initial,
        P_down_initial=P_initial,
        N=N,
        t_end=100.0,
    )

    if not np.allclose(
        result["P_final"],
        P_initial,
        rtol=0.0,
        atol=1.0e-8,
    ):
        raise AssertionError(
            "Uniform-pressure test failed: "
            "the final pressure field changed."
        )

    if not np.allclose(
        result["P_up"],
        P_initial,
        rtol=0.0,
        atol=1.0e-8,
    ):
        raise AssertionError(
            "Uniform-pressure test failed: "
            "upstream pressure changed."
        )

    if not np.allclose(
        result["P_down"],
        P_initial,
        rtol=0.0,
        atol=1.0e-8,
    ):
        raise AssertionError(
            "Uniform-pressure test failed: "
            "downstream pressure changed."
        )

    print("PASS: uniform-pressure test")


# --------------------------------------------------------------------------------------
# Test 2: Zero permeability
# --------------------------------------------------------------------------------------

def test_zero_permeability():
    """
    If intrinsic permeability is zero, no gas can move through the core.
    The timestep stability formula is bypassed and every pressure should
    remain at its initial value.
    """

    P_up_initial = 1.2e6
    P_down_initial = 1.0e6

    result = run_fdm(
        V_up=V_UP,
        V_down=V_DOWN,
        L=L,
        D=D,
        porosity=POROSITY,
        k=0.0,
        mu=MU,
        b=B,
        P_up_initial=P_up_initial,
        P_down_initial=P_down_initial,
        N=N,
        t_end=100.0,
    )

    expected_final = np.full(N + 1, P_down_initial)
    expected_final[0] = P_up_initial

    if not np.array_equal(
        result["P_final"],
        expected_final,
    ):
        raise AssertionError(
            "Zero-permeability test failed: "
            "pressure field changed."
        )

    if result["num_steps"] != 0:
        raise AssertionError(
            "Zero-permeability test failed: "
            "solver performed timesteps."
        )

    if result["dt"] is not None:
        raise AssertionError(
            "Zero-permeability test failed: "
            "dt should be None."
        )

    if result["max_inventory_error"] != 0.0:
        raise AssertionError(
            "Zero-permeability test failed: "
            "inventory changed."
        )

    print("PASS: zero-permeability test")


# --------------------------------------------------------------------------------------
# Test 3: Baseline regression
# --------------------------------------------------------------------------------------

def test_baseline_regression():
    """
    Verify that the reusable solver reproduces the previously verified
    N = 40 baseline solution after refactoring.
    """

    result = run_fdm(
        V_up=V_UP,
        V_down=V_DOWN,
        L=L,
        D=D,
        porosity=POROSITY,
        k=K,
        mu=MU,
        b=B,
        P_up_initial=1.2e6,
        P_down_initial=1.0e6,
        N=40,
        t_end=5000.0,
        safety_factor=0.8,
    )

    expected_P_up = 1.080270e6
    expected_P_down = 1.080269e6
    expected_inventory_error = 2.012071e-5

    if not np.isclose(
        result["P_up"][-1],
        expected_P_up,
        rtol=0.0,
        atol=1.0,
    ):
        raise AssertionError(
            "Baseline regression failed: "
            "final upstream pressure changed."
        )

    if not np.isclose(
        result["P_down"][-1],
        expected_P_down,
        rtol=0.0,
        atol=1.0,
    ):
        raise AssertionError(
            "Baseline regression failed: "
            "final downstream pressure changed."
        )

    if not np.isclose(
        result["max_inventory_error"],
        expected_inventory_error,
        rtol=1.0e-4,
        atol=0.0,
    ):
        raise AssertionError(
            "Baseline regression failed: "
            "inventory error changed."
        )

    print("PASS: baseline regression test")


# --------------------------------------------------------------------------------------
# Run verification tests
# --------------------------------------------------------------------------------------

if __name__ == "__main__":

    print()
    print("Running FDM verification tests...")
    print()

    test_uniform_pressure()
    test_zero_permeability()
    test_baseline_regression()

    print()
    print("All FDM verification tests passed.")