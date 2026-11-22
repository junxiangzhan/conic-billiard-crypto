# Security policy

This project is a **research prototype**. The maintainers do not recommend it for production cryptography, storing secrets, authentication, or key exchange.

In particular, modular square roots and affine group operations are variable-time; exception handling is incomplete for some projective points; no formal security analysis, protocol authentication, KDF, or side-channel testing has been performed. Some samples intentionally use trivially small groups.

Please avoid opening public issues containing secret material. For research discussion, provide only public parameters and reproducible non-secret test vectors.
