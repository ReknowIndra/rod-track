# rod-track
Tracking of rodlike particles in microscopic video sequences  (Python)

The algorithm uses a grid-search-based approach to find a predefined section in consecutive images. The originally selected region is translated and rotated in the vicinity of its previous position. As the object to track may change its appearance over time, one version is implemented where the search region is updated after each step. This approach has the disadvantage of error propagation in the form of dead reckoning. The second version looks for the original region in all subsequent steps. In the application here this proved to be the more appropriate approach. The Python code is kept simple without more attention on adaptability rather than performance.
