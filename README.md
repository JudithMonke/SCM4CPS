Detection Analysis of FDIA using χ2-detector and Euclidian distance based Detector:
The Evolution of the communication network in the Smart Grid represented
by the introduction of advanced metering infrastructure, wireless
sensors and cloud computing introduces additional harmful vulnerabilities
that a skilled attacker can exploit to perform targeted FDIAs. We
performed simulations to compare the χ2-detector and the Euclidean Distance
based detector proposed in [2] in detecting random and targeted
False Data Injection Attacks.

The authors [2] proved that a well implemented FDIA, according
to the definitions in paper [3], is not recognised by the most widely used
detection method together with the kalman filter, the χ2 detector.
Instead, a method based on Euclidean distance is more sensitive towards signal
changes, and is also capable to detect False Data Injction Attacks.
Not only it is more effective, but it is also more responsive and faster to trigger
an alarm with respect to the χ2 detector in case of faults or Random attacks.
On the other hand, this reactivity of the Euclidian detector makes it much
more susceptible to false alarms due to faults or noises, so in these cases the χ2
detector is preferable since it handles soft errors better [2].
Moreover, since the Euclidean detector reconstructs the signal from the state
estimates and compares it with the measured signal, is more resource intensive
than the χ2-detector, which instead only computes the residue vector.
