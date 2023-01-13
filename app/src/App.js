
import { React, useState } from "react";
import { Navbar, Container, Row, Col, Table, Form, Button, Spinner } from "react-bootstrap";
import { MapContainer, TileLayer, Marker, Popup } from "react-leaflet";
import MarkerClusterGroup from "@changey/react-leaflet-markercluster";
import "bootstrap/dist/css/bootstrap.min.css";
import "leaflet/dist/leaflet.css";
import iconMarker from "leaflet/dist/images/marker-icon.png";
import iconRetina from "leaflet/dist/images/marker-icon-2x.png";
import iconShadow from "leaflet/dist/images/marker-shadow.png";
import "@changey/react-leaflet-markercluster/dist/styles.min.css";
import L from "leaflet";

L.Icon.Default.mergeOptions({
  iconRetinaUrl: iconRetina,
  iconUrl: iconMarker,
  shadowUrl: iconShadow,
});

function App() {

  const [occurrences, setOccurrences] = useState([]);
  const [loading, setLoading] = useState(false);
  const [sequence, setSequence] = useState(`TAGTCATATGCTTGTCTCAAAGATAAGCCATGCATGTCTAAGTATAAGCGACTATACTGTGAAACTGCGA
ATGGCTCATTAAATCAGTTATGGTTTATTTGATGGTACCTTGCTACTTGGATAACCGTAGTAATTCTAGA
GCTAATACATGCAGGAGTTCCCGACTCACGGAGGGATGTATTTATTAGATAAGAAACCAAACCGGTCTCC
GGTTGCGTGCTGAGTCATAATAACTGCTCGAATCGCACGGCTCTACGCCGGCGATGGTTCATTCAAATTT
CTGCCCTATCAGCTTTCGATGGTAGGATAGAGGCCTACCATGGCGTTAACGGGTAACGGAGAATTAGGGT
TCGATTCCGGAGAGGGAGCCTGAGAAATGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTGCCC
GAATCCTGACACAGGGAGGTAGTGACAAGAAATAACAATACAGGGCTATTTTAGTCTTGTAATTGGAATG
AGTACAATTTACATCTCTTCACGAGGATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCC
AGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAACGCTCGTAGTCGGATTTCGGGGCGGGCCGACC
GGTCTGCCGATGGGTATGCACTGGCCGGCGCGTCCTTCCACCCGGAGACCGCGCCTACTCTTAACTGAGC
GGGCGCGGGAGACGGGTCTTTTACTTTGAAAAAATCAGAGTGTTTCAAGCAGGCAGTCGCTCTTGCATGG
ATTAGCATGGGATAATGAAATAGGACTCTGGTGCTATTTTGTTGGTTTCGAACACCGGAGTAATGATTAA
CAGGGACAGTCAGGGGCACTCGTATTCCGCCGAGAGAGGTGAAATTCTCAGACCAGCGGAAGACGAACCA
CTGCGAAAGCATTTGCCAGGGATGTTTTCACTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGA
TACCGTCGTAGTCTTAACCATAAACCATGCCGACTAGGGATTGGAGGATGTTCCATTTGTGACTCCTTCA
GCACCTTTCGGGAAACTAAAGTCTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAA
TTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGG
TCCAGCACATTGTGAGGATTGACAGATTGAGAGCTCTTTCTTGATTCGATGGGTGGTGGTGCATGGCCGT
TCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCGCAGCCTGCTAAATAGCGA
CGCGAACCCTCCGTTCGCTGGAGCTTCTTAGAGGGACAACTTGTCTTCAACAAGTGGAAGTTCGCGGCAA
TAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGCACTCAACGAGTCTA
TCACCTTGACCGAGAGGTCCGGGTAATCTTTTGAAATTGCATCGTGATGGGGATAGATTATTGCAACTAT
TAATCTTCAACGAGGAATTCCTAGTAAGCGTGTGTCATCAGCGCACGTTGATTACGTCCCTGCCCTTTGT
ACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAGGCCCCCGGACTGCGGCGCCGCAGCTGGT
TCTCCAGCCGCGACGCCGCGGGAAGCTGTCCGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGG`);

  function handleSequenceChange(e) {
    setSequence(e.target.value);
  }

  function handleSearch() {
    setLoading(true);
    async function search() {
      const res = await fetch("https://api.sequence.obis.org/search?sequence=" + sequence.trim());
      const data = await res.json();
      setOccurrences(data["results"]);
      setLoading(false);
    }
    search(sequence);
  }

  return (
    <div className="App">
      <Navbar bg="light" expand="lg">
        <Container>
          <Navbar.Brand href="/">
            OBIS sequence search
          </Navbar.Brand>
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
            occurrences.map(occ => <Marker key={occ.occurrence_id} position={[occ.decimalLatitude, occ.decimalLongitude]} >
              <Popup>
                <Table className="text-sm table-sm">
                  <tbody>
                    <tr><td>ID</td><td><a href={"https://obis.org/occurrence/" + occ.occurrence_id} rel="noreferrer" target="_blank">{occ.occurrence_id}</a></td></tr>
                    <tr><td>dataset ID</td><td><a href={"https://obis.org/dataset/" + occ.dataset_id} rel="noreferrer" target="_blank">{occ.dataset_id}</a></td></tr>
                    <tr><td>scientificName</td><td>{occ.scientificName}</td></tr>
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

      </MapContainer>

      <Container className="mt-4 mb-3">
        <Row className="mt-3 mb-3">
          <Col lg={true} className="mb-3">
            <Form.Group className="mb-3" controlId="sequence">
              <Form.Label>Sequence</Form.Label>
              <Form.Control className="font-monospace" as="textarea" rows={8} value={sequence} onChange={handleSequenceChange} />
            </Form.Group>
            <Button variant="primary" onClick={handleSearch}>Search
            { loading &&
              <Spinner className="mx-2" animation="border" size="sm"></Spinner>
            }
            </Button>
            
          </Col>
        </Row>
        <Row>
          <Col>
            {
              occurrences.length ?
              <Table className="text-sm table-sm">
                <thead>
                  <tr>
                    <th>scientificName</th>
                    <th>alignment score</th>
                    <th>phylum</th>
                    <th>class</th>
                    <th>order</th>
                    <th>family</th>
                    <th>genus</th>
                    <th>dataset</th>
                  </tr>
                </thead>
                <tbody>
                  { occurrences.map((occ) => <tr key={occ.id}>
                    <td>{occ.scientificName}</td>
                    <td>{occ.as}</td>
                    <td>{occ.phylum}</td>
                    <td>{occ.class}</td>
                    <td>{occ.order}</td>
                    <td>{occ.family}</td>
                    <td>{occ.genus}</td>
                    <td><a href={"https://obis.org/dataset/" + occ.dataset_id} rel="noreferrer" target="_blank">link</a></td>
                  </tr>) }
                </tbody>
              </Table> :
              <p></p>
            }
          </Col>
        </Row>
      </Container>
    </div>
  );
}

export default App;
