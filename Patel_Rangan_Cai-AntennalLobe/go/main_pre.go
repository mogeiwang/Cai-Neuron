package main
import (
    "fmt"
    "math"
)

type (
    Integer int64
    Real float64
    Complex complex128
    Boolean bool
    String string
    Character rune
)

const (
// parameters for leak current
    C_m_PN Real = 1.0
    g_L_PN Real = 0.3
    E_L_PN Real = -64
    C_m_LN Real = 1.0
    g_L_LN Real = 0.3
    E_L_LN Real = -50

// parameters for I_Intrinsic
    g_K_PN   Real  = 3.6
    g_K_LN   Real  = 36
    g_Na_PN  Real  = 120
    g_A_PN   Real  = 1.43
    g_Ca_LN  Real  = 5.0
    g_CaK_LN Real  = 0.045

    E_K_PN   Real  = -87
    E_K_LN   Real  = -95
    E_Na_PN  Real  = 40
    E_A_PN   Real  = E_K_PN
    E_Ca_LN  Real  = 140
    E_CaK_LN Real  = E_K_LN
)

type xNdata struct { //PN and LN, together
    ID Integer
    genre Integer // LN (0), or PN (1)?
    t Real // current time
    x Real // position in space
    y Real
    z Real

    V Real // V_PN or V_LN
    C_m Real
    g_L Real
    E_L Real

    g_K Real // PN, LN
    g_Na Real //PN
    g_A Real // PN
    g_Ca Real // LN
    g_CaK Real //LN

    E_K Real // PN, LN
    E_Na Real //PN
    E_A Real // PN
    E_Ca Real // LN
    E_CaK Real //LN

    I_K Real // PN, LN
    I_Na Real //PN
    I_A Real // PN
    I_Ca Real // LN
    I_CaK Real //LN

    I_GABA Real
    I_nACH Real
    I_stim Real
    I_slow Real //PN
}

type xNfunc interface { // just collection, not used
    init(i,g Integer) Boolean
    update_I_K()    Real //PN, LN
    update_I_Na()   Real //PN
    update_I_A()    Real //PN
    update_I_Ca()   Real //LN
    update_I_CaK()  Real //LN

    update_I_GABA() Real
    update_I_nACH() Real
    update_I_stim() Real
    update_I_slow() Real //PN
}

func (x *xNdata) init(i,g Integer) Boolean {
    if g == 0 {//LN
        x.ID = i
        x.genre = g
        x.C_m = C_m_LN
        x.g_L = g_L_LN
        x.E_L = E_L_LN

        x.g_K   = g_K_LN
        x.g_Ca  = g_Ca_LN
        x.g_CaK = g_CaK_LN
        x.E_K   = E_K_LN
        x.E_Ca  = E_Ca_LN
        x.E_CaK = E_CaK_LN

    } else if g == 1 {//PN
        x.ID = i
        x.genre = g
        x.C_m = C_m_PN
        x.g_L = g_L_PN
        x.E_L = E_L_PN

        x.g_K   = g_K_PN
        x.g_Na  = g_Na_PN
        x.g_A   = g_A_PN
        x.E_K   = E_K_PN
        x.E_Na  = E_Na_PN
        x.E_A   = E_A_PN
    } else {
        return false
    }
    return true
}

func (x *xNdata) update_I_Intrinsic(s string) Boolean {
    if x.genre == 0 {//LN
        if s == "K" {
        } else if s == "Ca" {
        } else if s == "CaK" {
        }
    } else if x.genre == 1 {//PN
        if s == "K" {
        } else if s == "Na" {
        } else if s == "A"  {
        }
    } else {
        return false
    }
    return true
}

func main() {
    fmt.Print("hello, world \t", math.Pi)
}
