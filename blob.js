/** 
 * StarterCode for "Attack of the Blobs!" 
 * CS248B Fundamentals of Computer Graphics: Animation & Simulation
 * 
 * Fill in the the missing code (see TODO items).
 * Try reducing MAX_BLOBS to 1 to get started. 
 * Good luck!!
 * 
 * @author Doug L. James <djames@cs.stanford.edu> 
 * @date 10/28/2022
 */

const MAX_BLOBS = 100; /// TODO: 100 or more to complete "Attack of the Blobs!" challenge. Use just a few for testing. 
const DRAW_BLOB_PARTICLES = true;

const STIFFNESS_STRETCH = 1.0; // TODO: Set as you wish
const STIFFNESS_BEND = 1.0; //    TODO: Set as you wish
const STIFFNESS_AREA = 1.0; //    TODO: Set as you wish

const WIDTH = 1024;
const HEIGHT = 1024;
const PARTICLE_RADIUS = WIDTH / 400.0; // for rendering
const PARTICLE_MASS = 1.0;
const BLOB_PARTICLES = 15; // 12 (F22)
const BLOB_RADIUS = WIDTH / 25;

//////// IMPORTANT ARRAYS OF THINGS /////////
let particles = []; // All particles in the scene (rigid + blobs)
let edges = []; //     All edges in the scene (rigid + blobs)
let blobs = []; //     All blobs in the scene (increases over time)
let environment; //    Environment with all rigid edges available as getEdges()

let isPaused = true;
let nTimesteps = 0; // #frame-length timesteps taken
let detectedEdgeEdgeFailure = false; // Halts simulation and turns purple if true -- blobs win!

// Graph paper texture map:
let bgImage;

function preload() {
	bgImage = loadImage('graphpaper.jpg');
}

function setup() {
	createCanvas(WIDTH, HEIGHT);
	background(100);
	ellipseMode(RADIUS);
	environment = new Environment();
	//print("|particles|=" + particles.length + ",  |edge|=" + edges.length + ",  |blobs|=" + blobs.length);
}

/// Timesteps (w/ substeps) and draws everything.
function draw() {

	///// SIMULATE /////
	if (!isPaused) {
		// CREATE BLOBS 
		if (nTimesteps % 10 == 0) { 
			if (blobs.length < MAX_BLOBS) 
				createRandomBlob(); // tries to create one if free space available
		}

		// TIMESTEP!
		let dtFrame = 0.01;
		let nSubsteps = 1; // #times to split dtFrame
		for (let step = 0; step < nSubsteps; step++)
			advanceTime(dtFrame / nSubsteps);
		nTimesteps++;
	}

	///// RENDER /////
	push();
	background(0);
	environment.draw();
	for (let blob of blobs)
		blob.draw();
	pop();
	drawMouseForce();

	/// TEXT OUTPUT:
	push();
	textSize(18);
	noStroke();
	fill(0);
	text("#BLOBS: " + blobs.length, 10, 20);
	text("#EDGES: " + edges.length, 10, 40);
	text("#PARTICLES: " + particles.length, 10, 60);
	pop();
}

function keyPressed() {
	if (keyCode == 32) // spacebar
		isPaused = !isPaused;
}

function advanceTime(dt) {
	environment.advanceTime(dt);

	//////////////////////////////////////
	////// GATHER PARTICLE FORCES ////////
	{
		// Clear forces:
		for (let particle of particles)
			particle.f.set(0, 0);

		gatherParticleForces_Gravity();

		// Damping (springs or otherwise -- you can add some if you want): 

		// Blob springs: 
		for (let blob of blobs) {
			blob.gatherForces_Stretch();
			blob.gatherForces_Bend();
			blob.gatherForces_Area();
		}

		gatherParticleForces_Penalty();

		// Mouse force (modify if you want):
		applyMouseForce();
	}

	//////////////////////////////////////////
	// Update velocity (using mass filtering):
	for (let particle of particles)
		acc(particle.v, dt * particle.invMass(), particle.f)

	//////////////////////////////////////////
	// Collision filter: Correct velocities //
	applyPointEdgeCollisionFilter();

	//////////////////////////////////////////
	// Update positions:
	for (let particle of particles)
		acc(particle.p, dt, particle.v)
}

function applyPointEdgeCollisionFilter() {
	// TEMP HACK (remove!): rigid bounce off walls so they don't fly away
	for (let blob of blobs) blob.nonrigidBounceOnWalls();

	// TODO: Process all point-edge CCD impulses 
	// FIRST: Just rigid edges.
	let edgesToCheck = environment.getEdges();
	// SECOND: All rigid + blob edges (once you get this ^^ working)
	// edgesToCheck = edges;

	// Initially just brute force all-pairs checks, later use bounding volumes or better broad phase.

}

// Efficiently checks that no pair of edges overlap, where the pairs do not share a particle in common.
function verifyNoEdgeEdgeOverlap() {
	if (detectedEdgeEdgeFailure) return; // already done

	// TODO: Optional: Make faster with broad phase
	// SIMPLE: Brute force check on edges i<j:
	for (let i = 0; i < edges.length - 1; i++) {
		let ei = edges[i];
		for (let j = i + 1; j < edges.length; j++) {
			let ej = edges[j];
			if (checkEdgeEdgeOverlap(ei, ej)) {
				// HALT!
				detectedEdgeEdgeFailure = true;
				isPaused = true;
				return;
			}
		}
	}
}

function checkEdgeEdgeOverlap(ei, ej) {
	// TODO: Implement robust-enough test. 
	return false;
}

// Computes penalty forces between all point-edge pairs
function gatherParticleForces_Penalty() {

	let warmup = true;
	if (warmup) { // First just consider rigid environment edges:
		for (let edge of environment.getEdges()) {
			// TODO (part1): Apply point-edge force (if pt not on edge!)
		}
	} else { // Consider all rigid + blob edges:
		for (let edge of edges) {
			// TODO (part2): Apply point-edge force (if pt not on edge!)
		}
	}
}

function gatherParticleForces_Gravity() {
	let g = vec2(0, 100); //grav accel
	for (let particle of particles)
		acc(particle.f, particle.mass, g); // f += m g
}

// Blob currently being dragged by mouse forces, or undefined.
let mouseBlob;
// Selects closest blob for mouse forcing (mask out if using a GUI)
function mousePressed() {
	if (blobs.length == 0 || isPaused) return;

	// Set mouseBlob to blob closest to the mouse:
	let m = vec2(mouseX, mouseY);
	let minDist = 1000000000;
	let minCOM;
	let minBlob;
	for (let blob of blobs) {
		let com = blob.centerOfMass();
		if (com.dist(m) < minDist) {
			minDist = com.dist(m);
			minBlob = blob;
			minCOM = com;
		}
	}
	mouseBlob = minBlob;
}

function mouseReleased() {
	mouseBlob = undefined;
}

// Applies spring + damping force to all mouseBlob particles
function applyMouseForce() {
	if (mouseIsPressed && mouseBlob) {
		if (blobs.length < 1) return;
		let m = vec2(mouseX, mouseY);
		let blobCOM = mouseBlob.centerOfMass();
		let blobDist = blobCOM.dist(m);
		let mforce = sub(m, blobCOM).normalize().mult(100 * clamp(blobDist, 0, 100));


		// Apply force to blob particles:
		let P = mouseBlob.blobParticles();
		for (let part of P) {
			part.f.add(mforce);
			acc(part.f, -10.0, part.v); //some damping
		}
	}
}

// Draws line from the mouse to any forced mouseBlob
function drawMouseForce() {
	if (mouseIsPressed && mouseBlob) {
		if (blobs.length < 1) return;
		let m = vec2(mouseX, mouseY);
		let blobCOM = mouseBlob.centerOfMass();
		push();
		stroke(0);
		strokeWeight(5);
		line(m.x, m.y, blobCOM.x, blobCOM.y);
		pop();
	}
}

// Creates a default particle and adds it to particles list
function createParticle(x, y) {
	let p = new Particle(vec2(x, y), 1.0, PARTICLE_RADIUS);
	particles.push(p);
	return p;
}

class Particle {
	constructor(pRest, mass, radius) {
		this.pRest = vec2(pRest.x, pRest.y);
		this.p = vec2(pRest.x, pRest.y);
		this.v = vec2(0, 0);
		this.pin = false; // true if doesn't respond to forces
		this.mass = mass;
		this.radius = radius;
		this.f = vec2(0, 0);
	}
	invMass() {
		return (this.pin ? 0.0 : 1.0 / this.mass);
	}
	// Emits a circle
	draw() {
		// nobody = (this.pin ? fill("red") : fill(0)); // default colors (red if pinned)
		circle(this.p.x, this.p.y, this.radius); //ellipseMode(RADIUS);
	}
}

// Creates edge and adds to edge list
function createEdge(particle0, particle1) {
	let edge = new Edge(particle0, particle1);
	edges.push(edge);
	return edge;
}

// Edge spring
class Edge {
	// Creates edge spring of default stiffness, STIFFNESS_STRETCH
	constructor(particle0, particle1) {
		this.q = particle0;
		this.r = particle1;
		this.restLength = this.q.pRest.dist(this.r.pRest);
		this.stiffness = STIFFNESS_STRETCH;
	}
	// True if both particles are pinned
	isRigid() {
		return (this.q.pin && this.r.pin);
	}
	// Current length of edge spring
	length() {
		return this.q.p.dist(this.r.p);
	}
	// Rest length of edge spring
	lengthRest() {
		return this.restLength;
	}
	// Draws the unstylized line 
	draw() {
		let a = this.q.p;
		let b = this.r.p;
		line(a.x, a.y, b.x, b.y);
	}
}

// RIGID ENVIRONMENT COMPOSED OF LINE SEGMENTS (pinned Edges)
class Environment {

	constructor() {
		this.envParticles = [];
		this.envEdges = [];

		///// BOX /////
		let r = PARTICLE_RADIUS;
		this.p00 = createParticle(r, r);
		this.p01 = createParticle(r, HEIGHT - r);
		this.p11 = createParticle(WIDTH - r, HEIGHT - r);
		this.p10 = createParticle(WIDTH - r, r);
		this.p00.pin = this.p01.pin = this.p11.pin = this.p10.pin = true;
		this.envParticles.push(this.p00);
		this.envParticles.push(this.p01);
		this.envParticles.push(this.p11);
		this.envParticles.push(this.p10);
		this.envEdges.push(createEdge(this.p00, this.p01));
		this.envEdges.push(createEdge(this.p01, this.p11));
		this.envEdges.push(createEdge(this.p11, this.p10));
		this.envEdges.push(createEdge(this.p10, this.p00));

		///// OBSTACLES FOR FUN /////
		{
			// (F22) ANGLED LINES: 
			//for (let i = 0.5; i < 4; i++) this.createEnvEdge(i * width / 5, height / 2, (i + 1) * width / 5, height * 0.75);

			// (F23) PACHINKO PEGS:
			let n = 8;
			for (let i = 1; i < n; i++) {
				for (let j = 1; j < n; j++) {
					if ((i + j) % 2 == 0) { // alternating pegs
						//this.createPachinkoPeg(width * (i / n), height / 5 + height / 4 * 3 * (j / n), 11); // round 4-edge peg
						this.createPachinkoWedge(width * (i / n), height / 5 + height / 4 * 3 * (j / n), 11); // cheap 2-edge peg
					}
				}
			}
		}
	}

	// Returns all rigid-environment edges.
	getEdges() {
		return this.envEdges;
	}

	// Creates a lone rigid edge.
	createEnvEdge(x0, y0, x1, y1) {
		let p0 = createParticle(x0, y0);
		let p1 = createParticle(x1, y1);
		p0.pin = true;
		p1.pin = true;
		let e = createEdge(p0, p1);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envEdges.push(e);
	}
	// Create a lone roundish peg at (x0,y0) with radius r.
	createPachinkoPeg(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y + r);
		let p2 = createParticle(x - r, y);
		let p3 = createParticle(x, y - r);
		p0.pin = p1.pin = p2.pin = p3.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		let e23 = createEdge(p2, p3);
		let e30 = createEdge(p3, p0);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envParticles.push(p3);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
		this.envEdges.push(e23);
		this.envEdges.push(e30);
	}
	// Create a lighter-weight wedge-shaped peg at (x0,y0) with radius r.
	createPachinkoWedge(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y - r);
		let p2 = createParticle(x - r, y);
		p0.pin = p1.pin = p2.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
	}

	// Updates any moveable rigid elements
	advanceTime(dt) {}

	// Makes popcorn <jk> no it doesn't... 
	draw() {
		push();
		image(bgImage, 0, 0, WIDTH, HEIGHT);

		if (detectedEdgeEdgeFailure) { // HALT ON OVERLAP + DRAW PURPLE SCREEN
			push();
			fill(191, 64, 191, 150);
			rect(0, 0, width, height);
			pop();
		}

		stroke("blue");
		strokeWeight(PARTICLE_RADIUS);
		for (let edge of this.envEdges) {
			edge.draw();
		}
		fill("blue");
		noStroke();
		for (let particle of this.envParticles) {
			particle.draw();
		}
		pop(); // wait, it does pop :/ 
	}
}

// Creates a blob centered at (x,y), and adds things to lists (blobs, edges, particles).
function createBlob(x, y) {
	let b = new Blob(vec2(x, y));
	blobs.push(b);
	return b;
}

// Tries to create a new blob at the top of the screen. 
function createRandomBlob() {
	for (let attempt = 0; attempt < 5; attempt++) {
		let center = vec2(random(2 * BLOB_RADIUS, WIDTH - 2 * BLOB_RADIUS), BLOB_RADIUS * 1.3); //random horizontal spot
		// CHECK TO SEE IF NO BLOBS NEARBY:
		let tooClose = false;
		for (let blob of blobs) {
			let com = blob.centerOfMass();
			if (com.dist(center) < 3 * blob.radius) // too close
				tooClose = true;
		}
		// if we got here, then center is safe:
		if (!tooClose) {
			createBlob(center.x, center.y);
			return;
		}
	}
}

class Blob {
	constructor(centerRest) {
		this.radius = BLOB_RADIUS;
		this.centerRest = centerRest; // original location

		// CREATE PARTICLES:
		this.BP = []; //blob particles
		this.n = BLOB_PARTICLES;
		let v0 = vec2(random(-100, 100), random(200, 220));
		for (let i = 0; i < this.n; i++) {
			let xi = this.radius * cos(i / this.n * TWO_PI) + centerRest.x;
			let yi = this.radius * sin(i / this.n * TWO_PI) + centerRest.y;
			let particle = createParticle(xi, yi);
			particle.v.set(v0);
			this.BP.push(particle);
		}

		// CREATE EDGES FOR STRETCH SPRINGS + COLLISIONS:
		this.BE = []; // blob edges
		for (let i = 0; i < this.n; i++) {
			let p0 = this.BP[i];
			let p1 = this.BP[(i + 1) % this.n];
			this.BE.push(createEdge(p0, p1));
		}

		// SETUP APPEARANCE/FACIAL ELEMENTS: 
		// TODO
		let dc = 26;
		this.fillColor = color([221 + random(-dc, dc), 160 + random(-dc, dc), 221 + random(-dc, dc), 255]); // ("Plum"); // 221, 160, 221
	}

	blobParticles() {
		return this.BP;
	}

	// Loops over blob edges and accumulates stretch forces (Particle.f += ...)
	gatherForces_Stretch() {
		let k = STIFFNESS_STRETCH;
		for (let edge of this.BE) {
			// TODO

		}
	}
	// Loops over blob particles and accumulates bending forces (Particle.f += ...)
	gatherForces_Bend() {
		let k = STIFFNESS_BEND;
		for (let i = 0; i < this.n; i++) {
			// FIND particles before (p0) and after (p2) particle i (p1): 
			//let p0 = ????
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];

			// TODO

		}
	}
	// Loops over blob particles and gathers area compression forces (Particle.f += ...)
	gatherForces_Area() {
		let k = STIFFNESS_AREA;
		// TODO 

	}

	// Center of mass of all blob particles
	centerOfMass() {
		let com = vec2(0, 0);
		for (let particle of this.BP)
			acc(com, 1 / this.BP.length, particle.p); // assumes equal mass
		return com;
	}

	// Center of velocity of all blob particles
	centerOfVelocity() {
		let cov = vec2(0, 0);
		for (let particle of this.BP)
			acc(cov, 1 / this.BP.length, particle.v); // assumes equal mass
		return cov;
	}

	// Something simple to keep rigid blobs inside the box:
	rigidBounceOnWalls() {
		let pos = this.centerOfMass();
		let vel = this.centerOfVelocity();

		let R = BLOB_RADIUS + PARTICLE_RADIUS;

		// Boundary reflection (only if outside domain AND still moving outward):
		if ((pos.x < R && vel.x < 0) ||
			(pos.x > width - R && vel.x > 0)) {
			for (let particle of this.BP)
				particle.v.x *= -0.4;
		}
		if ((pos.y < R && vel.y < 0) ||
			(pos.y > height - R && vel.y > 0)) {
			for (let particle of this.BP)
				particle.v.y *= -0.4;
		}
	}

	// Something simple to keep nonrigid blob particles inside the box:
	nonrigidBounceOnWalls() {
		let R = PARTICLE_RADIUS;
		for (let particle of this.BP) {
			let pos = particle.p;
			let vel = particle.v;

			// Boundary reflection (only if outside domain AND still moving outward):
			if ((pos.x < R && vel.x < 0) ||
				(pos.x > width - R && vel.x > 0)) {
				vel.x *= -0.4;
			}
			if ((pos.y < R && vel.y < 0) ||
				(pos.y > height - R && vel.y > 0)) {
				vel.y *= -0.4;
			}
		}
	}

	draw() {
		push();
		strokeWeight(PARTICLE_RADIUS);
		stroke("DarkOrchid"); //BlueViolet");
		fill(this.fillColor); { // draw blob
			beginShape(TESS);
			for (let particle of this.BP)
				vertex(particle.p.x, particle.p.y);
			endShape(CLOSE);
		}

		if (DRAW_BLOB_PARTICLES) {
			fill("DarkOrchid");
			for (let particle of this.BP)
				circle(particle.p.x, particle.p.y, PARTICLE_RADIUS);
		}

		this.drawBlobFace();
		pop();
	}

	drawBlobFace() {
		push();
		// TODO: Draw your face here! :D 

		// CENTER OF MASS eyeball for now :/
		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		stroke(0);
		fill(255);
		circle(com.x, com.y, 5 * PARTICLE_RADIUS);
		fill(0);
		circle(com.x + 0.01 * cov.x + 3 * sin(nTimesteps / 3), com.y + 0.01 * cov.y + random(-1, 1), PARTICLE_RADIUS);
		pop();
	}

}


/////////////////////////////////////////////////////////////////
// Some convenient GLSL-like macros for p5.Vector calculations //
/////////////////////////////////////////////////////////////////
function length(v) {
	return v.mag();
}

function dot(x, y) {
	return x.dot(y);
}

function dot2(x) {
	return x.dot(x);
}

function vec2(a, b) {
	return createVector(a, b);
}

function vec3(a, b, c) {
	return createVector(a, b, c);
}

function sign(n) {
	return Math.sign(n);
}

function clamp(n, low, high) {
	return constrain(n, low, high);
}

function add(v, w) {
	return p5.Vector.add(v, w);
}

function sub(v, w) {
	return p5.Vector.sub(v, w);
}

function absv2(v) {
	return vec2(Math.abs(v.x), Math.abs(v.y));
}

function maxv2(v, n) {
	return vec2(Math.max(v.x, n), Math.max(v.y, n));
}

function minv2(v, n) {
	return vec2(Math.min(v.x, n), Math.min(v.y, n));
}

function vertexv2(p) {
	vertex(p.x, p.y);
}

// v += a*w
function acc(v, a, w) {
	v.x += a * w.x;
	v.y += a * w.y;
}

function rotateVec2(v, thetaRad) {
	const c = cos(thetaRad);
	const s = sin(thetaRad);
	return vec2(c * v.x - s * v.y, s * v.x + c * v.y);
}
