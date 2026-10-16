import csv
from pathlib import Path
import numpy as np
import pytest
from virial_portfolio.model import R_KJ_MOL_K, VirialFit, fit_virial, load_pressure_csv, design_matrix


def example_data():
    u=np.r_[np.linspace(.05,2.5,25),np.linspace(.02,.7,20)]
    t=np.r_[np.full(25,273.),np.full(20,298.)]
    # Formula is explicit and separate from the implementation's prediction routine.
    z=(-2500+500*u-30*u*u+1.5*u**3)/t+12-u
    return u,t,z


def test_known_coefficients_recovered():
    f=fit_virial(*example_data())
    np.testing.assert_allclose(f.a,[-2500,500,-30,1.5],rtol=0,atol=2e-8)
    np.testing.assert_allclose(f.b,[12,-1],rtol=0,atol=2e-10)
    assert f.rmse_log<1e-12 and f.rank==6


def test_row_order_does_not_change_result():
    u,t,z=example_data();order=np.random.default_rng(19).permutation(len(u))
    a,b=fit_virial(u,t,z),fit_virial(u[order],t[order],z[order])
    np.testing.assert_allclose(np.r_[a.a,a.b],np.r_[b.a,b.b],rtol=0,atol=1e-8)


def test_common_pressure_scale_changes_only_b0_and_not_heat():
    u,t,z=example_data();a=fit_virial(u,t,z);b=fit_virial(u,t,z+np.log(760))
    np.testing.assert_allclose(a.a,b.a,rtol=0,atol=1e-8)
    assert b.b[0]-a.b[0]==pytest.approx(np.log(760),abs=1e-10)
    assert a.q0==pytest.approx(b.q0,abs=1e-10)


def test_hand_computed_prediction_and_heat():
    f=VirialFit(np.array([-2400.,80.]),np.array([10.,.1]),0,0,None,0,0)
    # N=1, T=300: z=-2400/300+80/300+10+0.1.
    z=2.366666666666667
    assert float(f.log_p_over_n(1,300))==pytest.approx(z,abs=1e-14)
    assert float(f.pressure(1,300))==pytest.approx(np.exp(z),abs=1e-12)
    assert float(f.heat(1))==pytest.approx(R_KJ_MOL_K*2320,abs=1e-12)


def test_heat_matches_fixed_loading_temperature_derivative():
    f=fit_virial(*example_data());u=.3;t1,t2=285.,285.001
    numerical=-R_KJ_MOL_K*(np.log(float(f.pressure(u,t2)))-np.log(float(f.pressure(u,t1))))/(1/t2-1/t1)
    assert float(f.heat(u))==pytest.approx(numerical,abs=2e-8)


def test_analytic_henry_inverse_and_zero_loading_limit():
    f=VirialFit(np.array([0.]),np.array([2.]),0,0,None,0,0)
    assert float(f.pressure(0,300))==0
    assert f.inverse_loading(10,300,(0,4))==pytest.approx(10/np.exp(2),abs=1e-11)


def test_general_inverse_round_trip():
    f=fit_virial(*example_data())
    for u in [.1,.6,1.7,2.2]:
        assert f.inverse_loading(float(f.pressure(u,273)),273,(.05,2.5))==pytest.approx(u,abs=1e-10)


def test_nonmonotonic_inverse_is_rejected():
    f=VirialFit(np.array([0.]),np.array([0.,-2.]),0,0,None,0,0)
    with pytest.raises(ValueError,match='strictly increasing'):f.inverse_loading(.1,300,(.01,2))


def test_interior_turning_region_is_detected_even_when_endpoints_are_positive():
    # g(u)=1-6u+6u^2 is positive at u=0 and u=1, negative in between.
    f=VirialFit(np.array([0.]),np.array([0.,-6.,3.]),0,0,None,0,0)
    with pytest.raises(ValueError,match='strictly increasing'):f.inverse_loading(.01,300,(0,1))


@pytest.mark.parametrize('value',[0,-1,np.nan,np.inf])
def test_invalid_training_loadings(value):
    u,t,z=example_data();u[0]=value
    with pytest.raises(ValueError):fit_virial(u,t,z)


@pytest.mark.parametrize('value',[0,-273,np.nan,np.inf])
def test_invalid_temperatures(value):
    u,t,z=example_data();t[0]=value
    with pytest.raises(ValueError):fit_virial(u,t,z)


@pytest.mark.parametrize('value',[np.nan,np.inf,-np.inf])
def test_invalid_log_response(value):
    u,t,z=example_data();z[0]=value
    with pytest.raises(ValueError):fit_virial(u,t,z)


def test_single_temperature_cannot_identify_temperature_terms():
    u,t,z=example_data()
    with pytest.raises(ValueError,match='rank-deficient'):fit_virial(u,np.full_like(t,273),z)


def test_constant_loading_cannot_identify_polynomial_terms():
    u,t,z=example_data()
    with pytest.raises(ValueError,match='rank-deficient'):fit_virial(np.ones_like(u),t,z)


def test_too_few_rows():
    u,t,z=example_data()
    with pytest.raises(ValueError,match='insufficient'):fit_virial(u[:5],t[:5],z[:5])


@pytest.mark.parametrize('m,n',[(True,1),(-1,1),(1.5,1),(3,9)])
def test_invalid_model_orders(m,n):
    with pytest.raises(ValueError):fit_virial(*example_data(),m=m,n=n)


def test_shape_mismatch_and_empty_arrays():
    u,t,z=example_data()
    with pytest.raises(ValueError):fit_virial(u,t[:-1],z)
    with pytest.raises(ValueError):fit_virial(u,t,z[:-1])
    with pytest.raises(ValueError):fit_virial([],[],[])
    with pytest.raises(ValueError):fit_virial(u[:,None],t,z)


def test_constant_response_reports_undefined_r2():
    u,t,z=example_data();f=fit_virial(u,t,np.full_like(z,3))
    assert f.r2_log is None


def test_prediction_overflow_is_explicit():
    f=VirialFit(np.array([0.]),np.array([1000.]),0,0,None,0,0)
    with pytest.raises(ValueError,match='overflow'):f.pressure(1,300)


def test_inverse_target_outside_domain():
    f=fit_virial(*example_data())
    with pytest.raises(ValueError,match='outside'):f.inverse_loading(1e12,273,(.05,2.5))


def test_pressure_and_loading_unit_conversion(tmp_path):
    u,t,z=example_data();p=u*np.exp(z)
    def save(name,pfactor,nfactor):
        dest=tmp_path/name
        with dest.open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['temperature_K','uptake','pressure']);w.writerows(zip(t,u/nfactor,p/pfactor))
        return dest
    torr=save('torr.csv',1,1)
    bar=save('bar.csv',76000000/101325,1000)
    a=fit_virial(*load_pressure_csv(torr))
    b=fit_virial(*load_pressure_csv(bar,pressure_unit='bar',loading_unit='mol/g'))
    np.testing.assert_allclose(np.r_[a.a,a.b],np.r_[b.a,b.b],rtol=0,atol=1e-8)


@pytest.mark.parametrize('pressure',[0,-1,'nan','inf'])
def test_invalid_csv_pressures(tmp_path,pressure):
    f=tmp_path/'bad.csv';f.write_text(f'temperature_K,uptake,pressure\n273,1,{pressure}\n')
    with pytest.raises(ValueError):load_pressure_csv(f)


def test_csv_missing_columns_and_unknown_units(tmp_path):
    p=tmp_path/'bad.csv';p.write_text('x,y\n1,2\n')
    with pytest.raises(ValueError):load_pressure_csv(p)
    with pytest.raises(ValueError):load_pressure_csv(p,pressure_unit='psi')


def test_published_a0_arithmetic_rounds_to_reported_q0():
    q=-R_KJ_MOL_K*(-2466.79553)
    assert q==pytest.approx(20.510079220812507,abs=1e-12)
    assert round(q)==21
