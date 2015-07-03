# This file is part of QuTiP: Quantum Toolbox in Python.
#
#    Copyright (c) 2011 and later, Paul D. Nation and Robert J. Johansson.
#    All rights reserved.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions are
#    met:
#
#    1. Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#
#    2. Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#
#    3. Neither the name of the QuTiP: Quantum Toolbox in Python nor the names
#       of its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################


from numpy.testing import assert_, assert_equal, run_module_suite

from qutip import rand_unitary, qeye, rand_herm, rand_dm, rand_ket
from qutip.sparse import sp_eigs


class TestRand:
    """
    A test class for the built-in random quantum object generators.
    """

    def testRandUnitary(self):
        "random Unitary"

        U = [rand_unitary(5) for k in range(5)]
        for u in U:
            assert_(u * u.dag() == qeye(5))

    def testRandherm(self):
        "random hermitian"

        H = [rand_herm(5) for k in range(5)]
        for h in H:
            assert_equal(h.isherm, True)

    def testRanddm(self):
        "random density matrix"

        R = [rand_dm(5) for k in range(5)]
        for r in R:
            assert_(r.tr() - 1.0 < 1e-15)
            # verify all eigvals are >=0
            assert_(not any(sp_eigs(r.data, r.isherm, vecs=False)) < 0)
            # verify hermitian
            assert_(r.isherm)

    def testRandket(self):
        "random ket"
        P = [rand_ket(5) for k in range(5)]
        for p in P:
            assert_equal(p.type == 'ket', True)

if __name__ == "__main__":
    run_module_suite()
