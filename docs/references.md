# References

The computational core of this project is derived directly from the law of reflection; the genus-1, Poncelet and discrete-logarithm parts build on existing work. The links below are entry points for reading.

1. P. Corvaja, D. Masser, U. Zannier, **Finiteness theorems on elliptical billiards and a variant of the dynamical Mordell–Lang conjecture**, *Proceedings of the London Mathematical Society* (2023). [Article](https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/plms.12561). Section 2.2 describes the genus-1 state curve for a fixed caustic and the group translation explicitly.
2. H. Stachel, **The geometry of billiards in ellipses and their Poncelet grids**, *Journal of Geometry* (2022). [Article](https://link.springer.com/article/10.1007/s00022-021-00606-2). Geometry of confocal caustics for ellipses and hyperbolas.
3. **Algebraic billiards in the Fermat hyperbola**, *Advances in Mathematics* (2025). [Article](https://www.sciencedirect.com/science/article/pii/S0001870825005729). General algebraic billiards, including the integrable cases with quadratic boundaries.
4. **Ivory's theorem revisited**, *Journal of Integrable Systems* (2017). [Article](https://academic.oup.com/integrablesystems/article/2/1/xyx006/3939821). Confocal quadrics and the geometry of reflection.
5. **A note on cyclic groups, finite fields, and the discrete logarithm problem** (1992). [Article](https://doi.org/10.1007/BF01189025). The DLP on genus-0 curves and on some matrix groups can reduce to the DLP in a finite field.

### Poncelet's general theorem (background for the pencil construction)

6. Ph. Griffiths, J. Harris, **On Cayley's explicit solution to Poncelet's porism**, *L'Enseignement Mathématique* 24 (1978), 31–40. Explains the Poncelet map as a translation on an elliptic curve.
7. H. J. M. Bos, C. Kers, F. Oort, D. W. Raven, **Poncelet's closure theorem**, *Expositiones Mathematicae* 5 (1987), 289–364. History of the theorem and several proofs.
8. V. Dragović, M. Radnović, **Poncelet Porisms and Beyond: Integrable Billiards, Hyperelliptic Jacobians and Pencils of Quadrics**, Birkhäuser, Frontiers in Mathematics (2011). Monograph on Poncelet's theorem, billiards and integrable systems.
9. M. Levi, S. Tabachnikov, **The Poncelet grid and billiards in ellipses**, *American Mathematical Monthly* 114(10) (2007), 895–908. The Poncelet grid: intersections of extended sides of a billiard orbit lie on confocal conics.

10. J. H. Silverman, **The Arithmetic of Elliptic Curves**, 2nd ed., Graduate Texts in Mathematics 106, Springer (2009). Section X.1: with full 2-torsion, the 2-descent map $P\mapsto x(P)-e_i$ is a group homomorphism (Prop. X.1.4, used in section 6 of [elliptic_correspondence.md](elliptic_correspondence.md)); III.3.1 and III.4.8 are used for the Weierstrass isomorphism and the group law.

The parts of these classical results that this project needs are re-proved by computer algebra in its own coordinates ([proofs.md](proofs.md)) and checked numerically over finite fields. Built on top of them here are the curve → billiard realizability condition (L11–L12), the turning rules at exceptional states, and the billiard translation of security conditions ([security_dictionary.md](security_dictionary.md)).
