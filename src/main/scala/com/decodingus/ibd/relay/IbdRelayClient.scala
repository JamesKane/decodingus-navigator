package com.decodingus.ibd.relay

import com.decodingus.ibd.crypto.EncryptedPayload
import com.decodingus.util.Logger
import io.circe.parser.decode
import io.circe.syntax.*
import sttp.client3.*
import sttp.ws.WebSocket

import java.util.concurrent.atomic.{AtomicBoolean, AtomicInteger}
import scala.concurrent.{ExecutionContext, Future, Promise}
import scala.util.{Failure, Success, Try}

/**
 * WebSocket relay client for IBD comparison data exchange.
 *
 * Connects to the AppView WebSocket relay endpoint, which forwards
 * encrypted payloads between two Navigator instances. The relay cannot
 * read the content — it only routes by session ID.
 *
 * @param relayUrl  Base WebSocket URL (e.g., "wss://decoding-us.com/api/v1/ibd/relay")
 * @param sessionId Unique session ID for this comparison
 * @param authToken Authentication token for relay access
 */
class IbdRelayClient(
                      relayUrl: String,
                      sessionId: String,
                      authToken: String
                    ):
  private val log = Logger[IbdRelayClient]
  private val connected = AtomicBoolean(false)
  private val retryCount = AtomicInteger(0)
  private val maxRetries = 5

  @volatile private var messageHandler: EncryptedPayload => Unit = _ => ()
  @volatile private var errorHandler: Throwable => Unit = _ => ()
  @volatile private var currentWebSocket: Option[WebSocket[Future]] = None

  def isConnected: Boolean = connected.get()

  def onMessage(handler: EncryptedPayload => Unit): Unit =
    messageHandler = handler

  def onError(handler: Throwable => Unit): Unit =
    errorHandler = handler

  /**
   * Connect to the relay WebSocket.
   * Returns a Future that completes when the connection is established.
   */
  def connect()(implicit ec: ExecutionContext): Future[Unit] =
    val wsUrl = s"$relayUrl?session=$sessionId&token=$authToken"
    val backend = HttpClientFutureBackend()

    basicRequest
      .get(uri"$wsUrl")
      .response(asWebSocketAlways[Future, Unit] { ws =>
        connected.set(true)
        retryCount.set(0)
        currentWebSocket = Some(ws)
        log.info(s"Connected to IBD relay for session $sessionId")
        receiveLoop(ws)
      })
      .send(backend)
      .map(_ => ())
      .recoverWith { case e: Exception =>
        connected.set(false)
        val attempt = retryCount.incrementAndGet()
        if attempt <= maxRetries then
          val delayMs = math.min(1000L * math.pow(2, attempt - 1).toLong, 30000L)
          log.warn(s"Relay connection failed (attempt $attempt/$maxRetries), retrying in ${delayMs}ms: ${e.getMessage}")
          Future {
            Thread.sleep(delayMs)
          }.flatMap(_ => connect())
        else
          log.error(s"Relay connection failed after $maxRetries attempts: ${e.getMessage}")
          Future.failed(e)
      }

  /**
   * Send an encrypted payload through the relay.
   */
  def send(payload: EncryptedPayload)(implicit ec: ExecutionContext): Future[Unit] =
    currentWebSocket match
      case Some(ws) =>
        val json = payload.asJson.noSpaces
        ws.sendText(json).recover { case e: Exception =>
          log.error(s"Failed to send payload: ${e.getMessage}")
          errorHandler(e)
        }
      case None =>
        Future.failed(new IllegalStateException("Not connected to relay"))

  /**
   * Disconnect from the relay.
   */
  def disconnect()(implicit ec: ExecutionContext): Future[Unit] =
    connected.set(false)
    currentWebSocket match
      case Some(ws) =>
        currentWebSocket = None
        ws.close().recover { case _ => () }
      case None =>
        Future.successful(())

  private def receiveLoop(ws: WebSocket[Future])(implicit ec: ExecutionContext): Future[Unit] =
    if !connected.get() then return Future.successful(())

    ws.receiveText().flatMap { text =>
      decode[EncryptedPayload](text) match
        case Right(payload) =>
          Try(messageHandler(payload)) match
            case Failure(e) =>
              log.error(s"Error in message handler: ${e.getMessage}")
              errorHandler(e)
            case Success(_) => ()
        case Left(error) =>
          log.warn(s"Failed to decode relay message: ${error.getMessage}")

      receiveLoop(ws)
    }.recoverWith { case e: Exception =>
      connected.set(false)
      currentWebSocket = None
      if retryCount.get() < maxRetries then
        log.warn(s"Relay connection lost, attempting reconnect: ${e.getMessage}")
        connect()
      else
        log.error(s"Relay connection lost permanently: ${e.getMessage}")
        errorHandler(e)
        Future.successful(())
    }
