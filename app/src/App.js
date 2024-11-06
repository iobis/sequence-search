
import { React, useState } from "react";
import { Navbar, Container, Row, Col, Table, Form, Button, Spinner, Nav } from "react-bootstrap";
import { MapContainer, TileLayer, Marker, Popup } from "react-leaflet";
import MarkerClusterGroup from "@changey/react-leaflet-markercluster";
import "bootstrap/dist/css/bootstrap.min.css";
import "leaflet/dist/leaflet.css";
import iconMarker from "leaflet/dist/images/marker-icon.png";
import iconRetina from "leaflet/dist/images/marker-icon-2x.png";
import iconShadow from "leaflet/dist/images/marker-shadow.png";
import "@changey/react-leaflet-markercluster/dist/styles.min.css";
import L from "leaflet";
import Control from "react-leaflet-custom-control";
import ReactSlider from "react-slider";

L.Icon.Default.mergeOptions({
  iconRetinaUrl: iconRetina,
  iconUrl: iconMarker,
  shadowUrl: iconShadow,
});

const examples = [
  `CCGCGGTAAGACGGAGGATGCAAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGGTGGTTTAATAAGTCAACTGTTAAATCTTGAGGCTCAACTTCAAAATCGCAGTCGAAACTATTAGACTAGAGTATAGTAGAGGTAAAGGGAATTTCCAGTGGAGCGGTGAAATGCGTAGATATTGGAAAGAACACCGATGGCGAAAGCACTTTACTGGGCTATTACTAACACTCAGAGACGAAAGCTAGGGTAG`,
  `AGCCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGGATTTCGGAGCGGGCCAGTTGGTCGGCCGCAAGGTCTGTTACTGACTGGTCTGCTCTTCTTCGCAAAGACTGCGTGTGCTCTTGACTGAGTGTGCGTAGGATTTACGACGTTTACTTTGAAAAAATTAGAGTGTTCAAAGCAGGC`,
  `TGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCTTAACCTGCTAAATAGTTACATGTAATTCCGGTTATGTGGGCAACTTCTTAGAGGGACTTTGTGTGTATAACGCAAGGAAGTTTGAGGCAATAACA`,
  `AGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCTAAAATGCCTGAATATTCGTGCATGGAATAATAGAATAGGATGTCGATCCTATTTTGTTGGTTTTCGGGACTCGACATAATGATAAATAGGGACAGTCGGGGGCATTTGTATTCAAACGACAGAGGTGAAATTCTTGGACCGTTTGAAGACAAACTACTGCGAAAGCATTTG`
]

function App() {

  const [occurrences, setOccurrences] = useState([]);
  const [table, setTable] = useState([]);
  const [loading, setLoading] = useState(false);
  const [slider, setSlider] = useState(95);
  const [maxSlider, ] = useState(100);
  const [exampleIndex, setExampleIndex] = useState(0);
  const [sequence, setSequence] = useState(examples[0]);

  function handleSequenceChange(e) {
    setSequence(e.target.value);
  }
  
  function handleLoadExample(e) {
    if (exampleIndex === examples.length - 1) {
      setExampleIndex(0);
      setSequence(examples[0]);
    } else {
      setExampleIndex(exampleIndex + 1);
      setSequence(examples[exampleIndex + 1]);
    }
  }
  
  function handleSliderChange(value) {
  setSlider(value);
}

function handleSearch() {
    setLoading(true);
    async function search() {
      // const res = await fetch("https://api.sequence.obis.org/search?sequence=" + sequence.trim());
      const res = await fetch("http://localhost:8000/search?sequence=" + sequence.trim());
      const data = await res.json();
      setOccurrences(data["results"]);
      setTable(data["table"]);
      setLoading(false);
    }
    search(sequence);
  }

  return (
    <div className="App d-flex flex-column min-h-100">
      <Navbar bg="light" expand="lg">
        <Container>
          <Navbar.Brand href="/">
            OBIS sequence search
          </Navbar.Brand>
          <Navbar.Toggle aria-controls="basic-navbar-nav" />
          <Navbar.Collapse id="basic-navbar-nav">
            <Nav className="ms-auto">
              <Nav.Link href="https://github.com/iobis/sequence-search" rel="noreferrer" target="_blank">GitHub</Nav.Link>
            </Nav>
          </Navbar.Collapse>
        </Container>
      </Navbar>

      <MapContainer
        id="map"
        center={[10, 0]}
        zoom={1}
        scrollWheelZoom={true}
      >
        <TileLayer
          attribution='&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
          url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
        />

        <MarkerClusterGroup maxClusterRadius={10}>
          {
            occurrences.filter(x => x.decimallatitude && x.decimallongitude && x.as > slider).map((occ, i) => <Marker key={i} position={[occ.decimallatitude, occ.decimallongitude]} >
              <Popup>
                <Table className="text-sm table-sm">
                  <tbody>
                    <tr><td>dataset ID</td><td><a href={"https://obis.org/dataset/" + occ.dataset_id} rel="noreferrer" target="_blank">{occ.dataset_id}</a></td></tr>
                    <tr><td>scientificName</td><td>{occ.scientificname}</td></tr>
                    <tr><td>phylum</td><td>{occ.phylum}</td></tr>
                    <tr><td>class</td><td>{occ.class}</td></tr>
                    <tr><td>order</td><td>{occ.order}</td></tr>
                    <tr><td>family</td><td>{occ.family}</td></tr>
                    <tr><td>genus</td><td>{occ.genus}</td></tr>
                    <tr><td>alignment score</td><td>{occ.as}</td></tr>
                  </tbody>
                </Table>
              </Popup>
            </Marker>)
          }
        </MarkerClusterGroup>

        { slider !== null &&
          <Control prepend position="bottomleft">
            <span>Identity threshold: { slider }%</span>
            <ReactSlider className="slider" thumbClassName="thumb" trackClassName="track" value={slider} min={80} max={maxSlider} onAfterChange={(value, index) =>
              handleSliderChange(value)
            } />
          </Control>
        }

      </MapContainer>

      <Container className="mt-4 mb-4">
        <Row className="mt-3 mb-3">
          <Col lg={true} className="mb-3" sm={8}>
            <Form.Group className="mb-3" controlId="sequence">
              <Form.Label><h5>Sequence</h5></Form.Label>
              <Form.Control className="font-monospace" as="textarea" rows={8} value={sequence} onChange={handleSequenceChange} />
            </Form.Group>
            <Button variant="primary" onClick={handleSearch}>Search
            { loading &&
              <Spinner className="mx-2" animation="border" size="sm"></Spinner>
            }
            </Button>
            <Button className="ms-2" variant="light" onClick={handleLoadExample}>Load example</Button>
          </Col>
          <Col sm={4}>
            <h5>About</h5>
            <p>This application aligns sequences against all sequence records in the OBIS database using <a href="https://github.com/torognes/vsearch" target="_blank" rel="noreferrer">VSEARCH</a>. Up to 100 distinct sequences as well as the associated occurrence records are returned ordered by sequence similarity. In the table below, occurrences are grouped by sequence, taxonomy, and dataset.</p>
          </Col>
        </Row>
        <Row>
          <Col>
            {
              table.length ?
              <Table className="text-sm table-sm">
                <thead>
                  <tr>
                    <th>scientificName</th>
                    <th>identity</th>
                    <th>phylum</th>
                    <th>class</th>
                    <th>order</th>
                    <th>family</th>
                    <th>genus</th>
                    <th>dataset</th>
                  </tr>
                </thead>
                <tbody>
                  { table.map((occ, i) => <tr key={i}>
                    <td>{occ.scientificname}</td>
                    <td>{occ.as} <span className="bar" style={{width: Math.max(0, (occ.as - 95) * 10)}}></span></td>
                    <td>{occ.phylum}</td>
                    <td>{occ.class}</td>
                    <td>{occ.order}</td>
                    <td>{occ.family}</td>
                    <td>{occ.genus}</td>
                    <td><a href={"https://obis.org/dataset/" + occ.dataset_id} rel="noreferrer" target="_blank">{occ.dataset_id}</a></td>
                  </tr>) }
                </tbody>
              </Table> :
              <p></p>
            }
          </Col>
        </Row>
      </Container>

      <footer className="mt-auto">
        <div className="footer mt-auto p-4 bg-dark">
          <Container>
            <Row>
              <div className="col-md">
                <p>
                <img alt="" className="footer-logo mx-2" src="/NEW UNESCO_logo_hor_white.png" />
                <img alt="" className="footer-logo mx-2" src="/lounsbery.png" />
                </p>
              </div>
              <div className="col-md text-white opacity-50">
                <p>This project is sponsored by the <a className="text-white" href="https://www.rlounsbery.org/" rel="noreferrer" target="_blank">Richard Lounsbery Foundation</a>.</p>
              </div>
            </Row>
          </Container>
        </div>
      </footer>
    </div>
  );
}

export default App;
