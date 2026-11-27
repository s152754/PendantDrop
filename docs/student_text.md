# Camera uncertainty in pendant-drop measurements: a short explanation

## What are we trying to do?
We measure surface tension from a photo of a hanging drop. The drop's shape depends on its surface
tension, so if we know the shape precisely, in millimetres, we can work out the surface tension.
A camera, however, gives us pixels, not millimetres. To convert between the two we need a camera
model, and that model is never perfect.
**Why this matters:** any error in the conversion changes the drop shape we measure, and therefore
the surface tension we calculate. This work answers: how much?

## The camera model has two parts
The **intrinsic** parameters describe the camera itself: how strongly the lens magnifies (focal
length), where the centre of the image is, and how the lens bends straight lines (distortion).
The **extrinsic** parameters describe where the camera is and how it is tilted, relative to a small
printed black-and-white square (an ArUco marker) mounted next to the drop.
**Why two parts:** they are measured in different ways, so their uncertainty is determined
differently too.

## How uncertain are the intrinsics?
We photographed a checkerboard many times and let OpenCV calibrate the camera from the corners of the
squares. The detected corners are never exactly right; they are off by a fraction of a pixel. So we
repeated the calibration 500 times, each time shifting the corners slightly at random, by the amount
we measured they are typically off. The spread of the 500 results tells us how uncertain each
intrinsic parameter is.
**Why this way:** it is hard to calculate directly how corner errors affect a calibration, but easy
to simply try it many times and look at the spread. This is called a bootstrap.

## How uncertain are the extrinsics?
The same idea, with the marker. The computer finds the marker's four corners and refines each one by
fitting the edges of the black square very precisely. From how well those edge fits line up, we know
how uncertain each corner is. Then we shift the corners randomly 1,000 times, and each time
calculate where the camera must be. The spread of those 1,000 positions is the extrinsic
uncertainty.

## Two assumptions about the corners, and why we may make them
**Within one corner**, the error in the x direction and in the y direction could be linked. The code
allows for that, but here the link is practically zero, because the marker sits almost straight in
the image: a corner's x position comes mostly from one edge and its y position from the other.
**Between corners**, we treat the errors as unrelated. That is allowed for random errors, because
each corner is fitted from its own stretch of the edges; no pixel is used by two corners.
**What this does not cover:** errors that move all corners together, for example if blur or ink
spread makes the whole marker look slightly bigger. A shared shift of 0.1 pixel turns out to be
small, but since one pixel is about 8 micrometres on the marker, a printed marker that is a little
larger than 10 mm could matter. That is why measuring the printed marker is on the to-do list.

## How does this reach the drop?
We take the drop outline in pixels and convert it to millimetres thousands of times, each time with a
slightly different camera model drawn from the uncertainties above. The outlines do not land in
exactly the same place every time. How much they spread is the uncertainty the camera adds to the
drop shape.
**Why a virtual drop:** for now, the outline comes from the drop model instead of a real photo.
That is fine for this question, because the camera error depends on where the outline is in the
image, not on what the image looks like.

## Independent or joint, and why it matters
The simple approach treats lens errors and position errors as unrelated: pick a random lens error,
pick a random position error, combine. But they are related. Imagine the lens magnifies slightly
more than we think. Then, to make the marker appear the size it does in the photo, the camera must be
slightly farther away than we think. The two errors partly cancel. The joint approach captures this
by recalculating the camera position for every lens error.
**The result:** with the joint approach the camera adds about 6 µm of uncertainty per point on the
drop outline; treated as independent it would be about 28 µm, so the simple approach is roughly
4–5 times too pessimistic.

## What is still simplified?
- The drop outline is computed from the model, not yet detected in a real photo.
- The noise in detecting the drop edge is an assumed value (0.3 pixel, about 2.4 µm).
- We assume the drop hangs exactly in the plane of the marker. This matters a lot: because the camera
  is only about 40 mm away, a drop just 0.27 mm in front of or behind the marker causes as much error
  as the whole camera model.