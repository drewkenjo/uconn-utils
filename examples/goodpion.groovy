import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.clas.physics.LorentzVector
import uconn.utils.pid.stefan.PionCandidate
import uconn.utils.pid.stefan.PionCandidate.Cut

def banknames = ['RUN::config', 'REC::Particle','REC::Traj']

args.each{fname->
  def reader = new HipoReader()
  reader.open(fname)
  def event = new Event()
  def factory = reader.getSchemaFactory()
  def banks = banknames.collect{new Bank(factory.getSchema(it))}

  while(reader.hasNext()) {
    reader.nextEvent(event)
    banks.each{event.read(it)}

    if(banks.every()) {
      def (runb,partb,trajb) = banks

      // all particles that pass pi+ cuts from PionCandidate class
      def pips0 = (0..<partb.getRows()).collect{ipart->PionCandidate.getPionCandidate(ipart, partb, trajb)}.findAll{it.ispip()}
      pips0.each{candidate->
        println("good pi+ is found")
        def pip = candidate.getLorentzVector()
        println("pi+ momentum = "+pip.p())
      }

      //specify the list of pi+ cuts you want to apply
      def pips1 = (0..<partb.getRows()).collect{ipart->PionCandidate.getPionCandidate(ipart, partb, trajb)}.findAll{it.ispip(Cut.PION_PID, Cut.DELTA_VZ)}
      pips1.each{candidate->
        println("good pi+ is found")
        def pip = candidate.getLorentzVector()
        println("pi+ momentum = "+pip.p())
      }

    }
  }

  reader.close()
}

